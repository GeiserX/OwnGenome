//! `og dev schema SCHEMA.json DOC.json`: does a JSON document have the
//! shape its JSON Schema states? It replaces `tests/schema/validate.py`.
//!
//! The keywords are the ones `tests/schema/summary.schema.json` uses:
//! `type` (a name or a list), `enum`, `required`, `properties`,
//! `additionalProperties` (a schema, `true` or `false`), `items`,
//! `minItems`, `maxItems`, `minimum`, `allOf` and local `$ref`
//! (`#/definitions/…`); `$schema`, `title`, `description` and
//! `definitions` are allowed. The whole schema is compiled before the
//! document is read, so a rule the validator does not understand anywhere
//! in it (an unknown keyword or type name, a `$ref` that leads nowhere) is
//! exit 2, never a pass.
//!
//! Output: `INVALID <path>: <what>` per error, then `DOC: N schema
//! error(s)` (exit 1) or `DOC: valid against SCHEMA` (exit 0). A broken
//! schema prints `SCHEMA ERROR: …` on stderr, exit 2. A document that is
//! not JSON prints `ERROR: …`, exit 1.
//!
//! Where it differs from the Python on purpose:
//! - `additionalProperties: false` refuses a key the schema does not list
//!   (the Python skipped it).
//! - `enum` compares by JSON type: `true` is not `1` (Python's `True == 1`
//!   let it pass); `1` and `1.0` are both numbers and equal.
//! - `minimum` refuses NaN (`nan < 0` is false in Python, so it passed).
//! - an unknown keyword, an unknown type name or a `$ref` to nothing is
//!   exit 2 wherever it is in the schema, also under a key the document
//!   does not have (the Python raised `KeyError`, exit 1, or never looked).

use indexmap::IndexMap;
use pycompat::json::{self, Value};
use pycompat::value::py_repr;
use std::collections::HashMap;
use std::io::Write;

const USAGE: &str = "usage: og dev schema SCHEMA.json DOC.json\n";

const KNOWN: [&str; 15] = [
    "$schema",
    "title",
    "description",
    "definitions",
    "type",
    "enum",
    "required",
    "properties",
    "additionalProperties",
    "items",
    "minItems",
    "maxItems",
    "minimum",
    "allOf",
    "$ref",
];

/// A JSON Schema type name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Type {
    Object,
    Array,
    String,
    Integer,
    Number,
    Boolean,
    Null,
}

impl Type {
    fn parse(s: &str) -> Option<Type> {
        Some(match s {
            "object" => Type::Object,
            "array" => Type::Array,
            "string" => Type::String,
            "integer" => Type::Integer,
            "number" => Type::Number,
            "boolean" => Type::Boolean,
            "null" => Type::Null,
            _ => return None,
        })
    }

    /// As the Python's `TYPES`: a bool is never an integer or a number.
    fn holds(self, v: &Value) -> bool {
        match self {
            Type::Object => matches!(v, Value::Object(_)),
            Type::Array => matches!(v, Value::Array(_)),
            Type::String => matches!(v, Value::Str(_)),
            Type::Integer => matches!(v, Value::Int(_)),
            Type::Number => matches!(v, Value::Int(_) | Value::Float(_)),
            Type::Boolean => matches!(v, Value::Bool(_)),
            Type::Null => matches!(v, Value::Null),
        }
    }
}

/// What a schema says about keys it does not list.
#[derive(Debug)]
enum Additional {
    Any,
    Schema(Box<Node>),
    Forbidden,
}

/// One compiled schema.
#[derive(Debug)]
struct Node {
    reference: Option<String>,
    all_of: Vec<Node>,
    /// The type names as written, for the message, and parsed.
    types: Option<(Vec<String>, Vec<Type>)>,
    enumeration: Option<Vec<Value>>,
    minimum: Option<f64>,
    minimum_text: String,
    required: Vec<String>,
    properties: IndexMap<String, Node>,
    additional: Additional,
    min_items: Option<usize>,
    max_items: Option<usize>,
    items: Option<Box<Node>>,
}

/// A schema compiled whole: the root and every `$ref` target.
struct Schema {
    root: Node,
    refs: HashMap<String, Node>,
}

fn compile_schema(raw: &Value) -> Result<Schema, String> {
    let mut wanted = Vec::new();
    let root = compile(raw, "$", &mut wanted)?;
    let mut refs = HashMap::new();
    while let Some((r, at)) = wanted.pop() {
        if refs.contains_key(&r) {
            continue;
        }
        let target = resolve(raw, &r).map_err(|e| format!("{at}: {e}"))?;
        let node = compile(target, &r, &mut wanted)?;
        refs.insert(r, node);
    }
    Ok(Schema { root, refs })
}

/// The node a local `$ref` points at.
fn resolve<'a>(root: &'a Value, r: &str) -> Result<&'a Value, String> {
    let Some(path) = r.strip_prefix("#/") else {
        return Err(format!("only local $ref is supported: {r}"));
    };
    let mut node = root;
    for part in path.split('/') {
        node = match node {
            Value::Object(o) => o.get(part),
            _ => None,
        }
        .ok_or_else(|| format!("$ref {r} points at nothing"))?;
    }
    Ok(node)
}

fn compile(v: &Value, at: &str, wanted: &mut Vec<(String, String)>) -> Result<Node, String> {
    let Value::Object(o) = v else {
        return Err(format!("{at}: a schema must be an object"));
    };
    let mut unknown: Vec<&str> = o
        .keys()
        .map(String::as_str)
        .filter(|k| !KNOWN.contains(k))
        .collect();
    if !unknown.is_empty() {
        unknown.sort_unstable();
        let list: Vec<String> = unknown
            .iter()
            .map(|k| repr(&Value::Str((*k).to_owned())))
            .collect();
        return Err(format!(
            "{at}: unsupported keyword(s) [{}]",
            list.join(", ")
        ));
    }
    let count = |k: &str| -> Result<Option<usize>, String> {
        match o.get(k) {
            None => Ok(None),
            Some(Value::Int(n)) if *n >= 0 => Ok(Some(usize::try_from(*n).unwrap_or(usize::MAX))),
            Some(_) => Err(format!("{at}: {k} must be a whole number")),
        }
    };
    let reference = match o.get("$ref") {
        None => None,
        Some(Value::Str(r)) => {
            wanted.push((r.clone(), at.to_owned()));
            Some(r.clone())
        }
        Some(_) => return Err(format!("{at}: $ref must be a string")),
    };
    let all_of = match o.get("allOf") {
        None => Vec::new(),
        Some(Value::Array(a)) => a
            .iter()
            .enumerate()
            .map(|(i, s)| compile(s, &format!("{at}.allOf[{i}]"), wanted))
            .collect::<Result<_, _>>()?,
        Some(_) => return Err(format!("{at}: allOf must be a list")),
    };
    let types = match o.get("type") {
        None => None,
        Some(t) => {
            let names: Vec<String> = match t {
                Value::Str(s) => vec![s.clone()],
                Value::Array(a) => a
                    .iter()
                    .map(|x| match x {
                        Value::Str(s) => Ok(s.clone()),
                        _ => Err(format!("{at}: type names must be strings")),
                    })
                    .collect::<Result<_, _>>()?,
                _ => return Err(format!("{at}: type must be a name or a list")),
            };
            let parsed = names
                .iter()
                .map(|n| Type::parse(n).ok_or_else(|| format!("{at}: unknown type {n}")))
                .collect::<Result<_, _>>()?;
            Some((names, parsed))
        }
    };
    let enumeration = match o.get("enum") {
        None => None,
        Some(Value::Array(a)) => Some(a.clone()),
        Some(_) => return Err(format!("{at}: enum must be a list")),
    };
    let (minimum, minimum_text) = match o.get("minimum") {
        None => (None, String::new()),
        Some(Value::Int(n)) => (Some(*n as f64), n.to_string()),
        Some(Value::Float(x)) if !x.is_nan() => (Some(*x), pycompat::fmt::repr(*x)),
        Some(_) => return Err(format!("{at}: minimum must be a number")),
    };
    let required = match o.get("required") {
        None => Vec::new(),
        Some(Value::Array(a)) => a
            .iter()
            .map(|x| match x {
                Value::Str(s) => Ok(s.clone()),
                _ => Err(format!("{at}: required must list strings")),
            })
            .collect::<Result<_, _>>()?,
        Some(_) => return Err(format!("{at}: required must be a list")),
    };
    let mut properties = IndexMap::new();
    match o.get("properties") {
        None => {}
        Some(Value::Object(p)) => {
            for (k, s) in p {
                properties.insert(k.clone(), compile(s, &format!("{at}.{k}"), wanted)?);
            }
        }
        Some(_) => return Err(format!("{at}: properties must be an object")),
    }
    let additional = match o.get("additionalProperties") {
        None | Some(Value::Bool(true)) => Additional::Any,
        Some(Value::Bool(false)) => Additional::Forbidden,
        Some(s) => Additional::Schema(Box::new(compile(s, &format!("{at}.*"), wanted)?)),
    };
    let items = match o.get("items") {
        None => None,
        Some(s) => Some(Box::new(compile(s, &format!("{at}[]"), wanted)?)),
    };
    if let Some(Value::Object(defs)) = o.get("definitions") {
        for (k, s) in defs {
            compile(s, &format!("#/definitions/{k}"), wanted)?;
        }
    } else if o.contains_key("definitions") {
        return Err(format!("{at}: definitions must be an object"));
    }
    Ok(Node {
        reference,
        all_of,
        types,
        enumeration,
        minimum,
        minimum_text,
        required,
        properties,
        additional,
        min_items: count("minItems")?,
        max_items: count("maxItems")?,
        items,
    })
}

/// Python's `type(v).__name__`.
fn type_name(v: &Value) -> &'static str {
    match v {
        Value::Null => "NoneType",
        Value::Bool(_) => "bool",
        Value::Int(_) => "int",
        Value::Float(_) => "float",
        Value::Str(_) => "str",
        Value::Array(_) => "list",
        Value::Object(_) => "dict",
    }
}

/// Python's `repr(v)`; JSON text where that port refuses (a character
/// above U+00FF inside a list or dict).
fn repr(v: &Value) -> String {
    py_repr(v).unwrap_or_else(|_| {
        json::dumps(
            v,
            &json::Dump {
                ensure_ascii: false,
                ..json::Dump::default()
            },
        )
        .unwrap_or_default()
    })
}

/// `f"{x!r:.60}"`: the repr cut at 60 characters.
fn repr60(v: &Value) -> String {
    repr(v).chars().take(60).collect()
}

/// Equal by JSON type: a bool never equals a number, an int equals a float
/// of the same value, objects compare as sets of keys.
fn same(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Int(x), Value::Float(y)) | (Value::Float(y), Value::Int(x)) => *x as f64 == *y,
        (Value::Array(x), Value::Array(y)) => {
            x.len() == y.len() && x.iter().zip(y).all(|(p, q)| same(p, q))
        }
        (Value::Object(x), Value::Object(y)) => {
            x.len() == y.len() && x.iter().all(|(k, p)| y.get(k).is_some_and(|q| same(p, q)))
        }
        _ => a == b,
    }
}

fn as_number(v: &Value) -> Option<f64> {
    match v {
        Value::Int(n) => Some(*n as f64),
        Value::Float(x) => Some(*x),
        _ => None,
    }
}

impl Schema {
    fn check(&self, node: &Node, v: &Value, path: &str, errors: &mut Vec<String>) {
        if let Some(r) = &node.reference {
            self.check(&self.refs[r], v, path, errors);
        }
        for sub in &node.all_of {
            self.check(sub, v, path, errors);
        }
        if let Some((names, types)) = &node.types
            && !types.iter().any(|t| t.holds(v))
        {
            errors.push(format!(
                "{path}: expected {}, got {} {}",
                names.join("/"),
                type_name(v),
                repr60(v)
            ));
            return;
        }
        if let Some(e) = &node.enumeration
            && !e.iter().any(|x| same(x, v))
        {
            errors.push(format!(
                "{path}: {} is not one of {}",
                repr60(v),
                repr(&Value::Array(e.clone()))
            ));
        }
        if let (Some(min), Some(x)) = (node.minimum, as_number(v)) {
            if x.is_nan() {
                errors.push(format!(
                    "{path}: nan is not a number >= minimum {}",
                    node.minimum_text
                ));
            } else if x < min {
                let shown = match v {
                    Value::Int(n) => n.to_string(),
                    _ => pycompat::fmt::repr(x),
                };
                errors.push(format!("{path}: {shown} < minimum {}", node.minimum_text));
            }
        }
        if let Value::Object(o) = v {
            for key in &node.required {
                if !o.contains_key(key) {
                    errors.push(format!("{path}: missing required key '{key}'"));
                }
            }
            for (key, sub) in &node.properties {
                if let Some(x) = o.get(key) {
                    self.check(sub, x, &format!("{path}.{key}"), errors);
                }
            }
            for (key, x) in o {
                if node.properties.contains_key(key) {
                    continue;
                }
                match &node.additional {
                    Additional::Any => {}
                    Additional::Schema(sub) => self.check(sub, x, &format!("{path}.{key}"), errors),
                    Additional::Forbidden => errors.push(format!(
                        "{path}: key '{key}' is not allowed (additionalProperties: false)"
                    )),
                }
            }
        }
        if let Value::Array(a) = v {
            if let Some(n) = node.min_items
                && a.len() < n
            {
                errors.push(format!("{path}: {} items < minItems {n}", a.len()));
            }
            if let Some(n) = node.max_items
                && a.len() > n
            {
                errors.push(format!("{path}: {} items > maxItems {n}", a.len()));
            }
            if let Some(sub) = &node.items {
                for (i, x) in a.iter().enumerate() {
                    self.check(sub, x, &format!("{path}[{i}]"), errors);
                }
            }
        }
    }
}

fn load(path: &str) -> Result<Value, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    let text = String::from_utf8(bytes).map_err(|_| format!("{path}: not UTF-8"))?;
    // A leading byte-order mark is an error, as in json.load.
    json::loads(&text).map_err(|e| match e {
        json::JsonError::Decode(what, at) => format!("{path}: not JSON: {what} at byte {at}"),
        other => format!("{path}: not JSON: {other:?}"),
    })
}

/// The errors of `doc` against `schema`, or why the schema is broken.
pub fn validate(schema: &Value, doc: &Value) -> Result<Vec<String>, String> {
    let s = compile_schema(schema)?;
    let mut errors = Vec::new();
    s.check(&s.root, doc, "$", &mut errors);
    Ok(errors)
}

/// Runs `og dev schema ARGS`; returns the exit code.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let [schema_path, doc_path] = args else {
        let _ = write!(err, "{USAGE}");
        return 2;
    };
    let schema = match load(schema_path) {
        Ok(s) => s,
        Err(e) => {
            let _ = writeln!(err, "SCHEMA ERROR: {e}");
            return 2;
        }
    };
    let compiled = match compile_schema(&schema) {
        Ok(c) => c,
        Err(e) => {
            let _ = writeln!(err, "SCHEMA ERROR: {e}");
            return 2;
        }
    };
    let doc = match load(doc_path) {
        Ok(d) => d,
        Err(e) => {
            let _ = writeln!(out, "ERROR: {e}");
            return 1;
        }
    };
    let mut errors = Vec::new();
    compiled.check(&compiled.root, &doc, "$", &mut errors);
    for e in &errors {
        let _ = writeln!(out, "INVALID {e}");
    }
    if errors.is_empty() {
        let _ = writeln!(out, "{doc_path}: valid against {schema_path}");
        0
    } else {
        let _ = writeln!(out, "{doc_path}: {} schema error(s)", errors.len());
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(schema: &str, doc: &str) -> Result<Vec<String>, String> {
        validate(&json::loads(schema).unwrap(), &json::loads(doc).unwrap())
    }

    #[test]
    fn additional_properties_false_refuses_an_extra_key() {
        let e = run(
            r#"{"type":"object","properties":{"a":{"type":"integer"}},"additionalProperties":false}"#,
            r#"{"a":1,"b":2}"#,
        )
        .unwrap();
        assert_eq!(
            e,
            ["$: key 'b' is not allowed (additionalProperties: false)"]
        );
        assert!(
            run(r#"{"additionalProperties":true}"#, r#"{"b":2}"#)
                .unwrap()
                .is_empty()
        );
        let e = run(
            r#"{"additionalProperties":{"type":"string"}}"#,
            r#"{"b":2}"#,
        )
        .unwrap();
        assert_eq!(e, ["$.b: expected string, got int 2"]);
    }

    #[test]
    fn enum_keeps_bool_and_number_apart() {
        let e = run(r#"{"enum":[1]}"#, "true").unwrap();
        assert_eq!(e, ["$: True is not one of [1]"]);
        assert!(run(r#"{"enum":[1]}"#, "1.0").unwrap().is_empty());
        assert!(
            run(r#"{"enum":[{"a":[1,"x"]}]}"#, r#"{"a":[1.0,"x"]}"#)
                .unwrap()
                .is_empty()
        );
        assert_eq!(run(r#"{"enum":[true]}"#, "1").unwrap().len(), 1);
    }

    #[test]
    fn minimum_refuses_nan() {
        let e = run(r#"{"type":"number","minimum":0}"#, "NaN").unwrap();
        assert_eq!(e, ["$: nan is not a number >= minimum 0"]);
        let e = run(r#"{"minimum":0.5}"#, "[0, -1]").unwrap();
        assert!(e.is_empty(), "minimum applies to numbers only");
        let e = run(r#"{"minimum":0.5}"#, "-1").unwrap();
        assert_eq!(e, ["$: -1 < minimum 0.5"]);
    }

    #[test]
    fn a_broken_schema_is_an_error_wherever_it_is() {
        let bad = |s: &str| run(s, "{}").unwrap_err();
        assert!(bad(r##"{"$ref":"#/definitions/nope"}"##).contains("points at nothing"));
        assert!(bad(r#"{"type":"integr"}"#).contains("unknown type integr"));
        assert!(
            bad(r#"{"properties":{"never":{"format":"date"}}}"#).contains("unsupported keyword")
        );
        assert!(bad(r#"{"definitions":{"x":{"pattern":"a"}}}"#).contains("unsupported keyword"));
        assert!(bad(r#"{"$ref":"other.json#/x"}"#).contains("only local"));
    }

    #[test]
    fn messages_are_the_python_ones() {
        let s = r##"{"definitions":{"n":{"type":"integer"}},"type":"object","required":["a","z"],
            "properties":{"a":{"type":"array","minItems":3,"maxItems":1,"items":{"$ref":"#/definitions/n"}},
            "s":{"enum":["x","y"]}}}"##;
        let e = run(s, r#"{"a":[1,"two",false],"s":"it's"}"#).unwrap();
        assert_eq!(
            e,
            [
                "$: missing required key 'z'",
                "$.a: 3 items > maxItems 1",
                "$.a[1]: expected integer, got str 'two'",
                "$.a[2]: expected integer, got bool False",
                "$.s: \"it's\" is not one of ['x', 'y']",
            ]
        );
        let long = format!("\"{}\"", "a".repeat(100));
        let e = run(r#"{"type":"integer"}"#, &long).unwrap();
        assert_eq!(
            e[0],
            format!("$: expected integer, got str '{}", "a".repeat(59))
        );
    }
}
