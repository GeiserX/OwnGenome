//! The single reader of PharmCAT's flat and source-nested gene maps.
use pycompat::json::{self, Value};
use std::collections::BTreeSet;
use std::io::Write;

trait Get {
    fn get(&self, key: &str) -> Option<&Value>;
}
impl Get for Value {
    fn get(&self, key: &str) -> Option<&Value> {
        if let Value::Object(o) = self {
            o.get(key)
        } else {
            None
        }
    }
}

fn diplotypes(g: &Value) -> Result<&[Value], String> {
    let empty = &[][..];
    let pick = ["sourceDiplotypes", "recommendationDiplotypes"]
        .into_iter()
        .filter_map(|k| g.get(k))
        .find(|v| pycompat::value::truthy(v));
    match pick {
        None => Ok(empty),
        Some(Value::Array(a)) => Ok(a),
        _ => Err("diplotypes must be an array".into()),
    }
}

/// All diplotypes and the genes with a named allele. Invalid structure is
/// unassessed, just like invalid JSON, and cannot produce an empty success.
pub fn read(data: &Value) -> Result<(BTreeSet<String>, BTreeSet<String>), String> {
    let Value::Object(root) = data else {
        return Err("report must be an object".into());
    };
    let mut rows = BTreeSet::new();
    let mut called = BTreeSet::new();
    let mut gene = |name: &str, g: &Value| -> Result<(), String> {
        for d in diplotypes(g)? {
            if !matches!(d, Value::Object(_)) {
                return Err("diplotype must be an object".into());
            }
            let mut names = Vec::with_capacity(2);
            for a in ["allele1", "allele2"] {
                let n = match d.get(a) {
                    None | Some(Value::Null) => "",
                    Some(Value::Object(o)) => match o.get("name") {
                        None | Some(Value::Null) => "",
                        Some(Value::Str(s)) => s,
                        _ => return Err("allele name must be a string".into()),
                    },
                    _ => return Err("allele must be an object".into()),
                };
                names.push(n);
            }
            rows.insert(format!("{name}\t{}", names.join("/")));
            if names.iter().any(|n| {
                !n.is_empty() && !["unknown", "none", "?"].contains(&n.to_lowercase().as_str())
            }) {
                called.insert(name.to_owned());
            }
        }
        Ok(())
    };
    match root.get("genes") {
        None | Some(Value::Null) => {}
        Some(Value::Object(genes)) => {
            for (key, val) in genes {
                let Value::Object(map) = val else { continue };
                if map.contains_key("sourceDiplotypes")
                    || map.contains_key("recommendationDiplotypes")
                {
                    gene(key, val)?;
                } else {
                    for (name, g) in map {
                        if matches!(g, Value::Object(_)) {
                            gene(name, g)?;
                        }
                    }
                }
            }
        }
        _ => return Err("genes must be an object".into()),
    }
    Ok((rows, called))
}

/// `og dev pharmcat diplotypes REPORT [--called|--called-lines]`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let valid = args.first().is_some_and(|s| s == "diplotypes")
        && (args.len() == 2
            || (args.len() == 3 && ["--called", "--called-lines"].contains(&args[2].as_str())));
    if !valid {
        let _ = writeln!(
            err,
            "usage: og dev pharmcat diplotypes REPORT [--called|--called-lines]"
        );
        return 2;
    }
    let result = std::fs::read_to_string(&args[1])
        .map_err(|e| e.to_string())
        .and_then(|s| json::loads(&s).map_err(|e| e.to_string()))
        .and_then(|v| read(&v));
    let (rows, called) = match result {
        Ok(r) => r,
        Err(e) => {
            let _ = writeln!(err, "unreadable report: {e}");
            return 0;
        }
    };
    let result = match args.get(2).map(String::as_str) {
        Some("--called") => writeln!(
            err,
            "{}",
            called.iter().cloned().collect::<Vec<_>>().join(" ")
        )
        .and_then(|()| writeln!(out, "{}", called.len())),
        Some("--called-lines") => {
            writeln!(out, "{}", called.into_iter().collect::<Vec<_>>().join("\n"))
        }
        _ => writeln!(out, "{}", rows.into_iter().collect::<Vec<_>>().join("\n")),
    };
    if result.is_ok() { 0 } else { 2 }
}
