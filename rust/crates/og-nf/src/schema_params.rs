//! The params of `nextflow_schema.json`: every property of every group
//! under `definitions` (or `$defs`), with its `default` when it has one.

use pycompat::json::Value;

/// `(name, default)` of each param, in file order. A group or a
/// `properties` that is not an object holds no params.
pub fn read(schema: &Value) -> Vec<(String, Option<Value>)> {
    let Value::Object(top) = schema else {
        return Vec::new();
    };
    let truthy = |k: &str| top.get(k).filter(|v| pycompat::value::truthy(v));
    let Some(Value::Object(groups)) = truthy("definitions").or_else(|| truthy("$defs")) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for group in groups.values() {
        let Value::Object(g) = group else { continue };
        let Some(Value::Object(props)) = g.get("properties") else {
            continue;
        };
        for (name, p) in props {
            let default = match p {
                Value::Object(p) => p.get("default").cloned(),
                _ => None,
            };
            out.push((name.clone(), default));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_params_and_defaults() {
        let s = pycompat::json::loads(
            r#"{"definitions": {}, "$defs": {"g": {"properties": {"a": {"default": 1}, "b": {}}}, "h": 3}}"#,
        )
        .unwrap();
        let p = read(&s);
        assert_eq!(p, [("a".into(), Some(Value::Int(1))), ("b".into(), None)]);
        assert!(read(&Value::Null).is_empty());
    }
}
