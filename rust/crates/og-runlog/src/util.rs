//! Small shared pieces: Python's path arithmetic, JSON building, and
//! Python's `==` on JSON values.

use indexmap::IndexMap;
use pycompat::json::Value;
use std::path::{Component, Path, PathBuf};

/// `os.path.normpath`: `.` and `a/..` removed lexically, no file system.
pub fn normpath(p: &Path) -> PathBuf {
    let mut out: Vec<Component<'_>> = Vec::new();
    for c in p.components() {
        match c {
            Component::CurDir => {}
            Component::ParentDir => match out.last() {
                Some(Component::Normal(_)) => {
                    out.pop();
                }
                Some(Component::RootDir) => {}
                _ => out.push(c),
            },
            _ => out.push(c),
        }
    }
    if out.is_empty() {
        return PathBuf::from(".");
    }
    out.iter().collect()
}

/// `os.path.abspath`: joined to the working folder, then [`normpath`].
pub fn abspath(p: &Path) -> PathBuf {
    if p.is_absolute() {
        normpath(p)
    } else {
        let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("/"));
        normpath(&cwd.join(p))
    }
}

/// `os.path.basename`.
pub fn basename(p: &Path) -> String {
    let s = p.to_string_lossy();
    s.rsplit('/').next().unwrap_or("").to_owned()
}

/// An object in insertion order.
pub type Obj = IndexMap<String, Value>;

pub fn s(v: impl Into<String>) -> Value {
    Value::Str(v.into())
}

pub fn opt_s(v: Option<impl Into<String>>) -> Value {
    v.map_or(Value::Null, |v| Value::Str(v.into()))
}

pub fn opt_f(v: Option<f64>) -> Value {
    v.map_or(Value::Null, Value::Float)
}

pub fn list(v: &[String]) -> Value {
    Value::Array(v.iter().map(|x| s(x.as_str())).collect())
}

/// Python's `a == b` for two values `json.load` made: numbers by value
/// (`1 == 1.0 == True`), dicts regardless of key order.
pub fn py_eq(a: &Value, b: &Value) -> bool {
    fn num(v: &Value) -> Option<f64> {
        match v {
            Value::Bool(b) => Some(f64::from(u8::from(*b))),
            Value::Int(i) => Some(*i as f64),
            Value::Float(x) => Some(*x),
            _ => None,
        }
    }
    match (a, b) {
        (Value::Null, Value::Null) => true,
        (Value::Int(x), Value::Int(y)) => x == y,
        (Value::Str(x), Value::Str(y)) => x == y,
        (Value::Array(x), Value::Array(y)) => {
            x.len() == y.len() && x.iter().zip(y).all(|(p, q)| py_eq(p, q))
        }
        (Value::Object(x), Value::Object(y)) => {
            x.len() == y.len() && x.iter().all(|(k, v)| y.get(k).is_some_and(|w| py_eq(v, w)))
        }
        _ => match (num(a), num(b)) {
            (Some(x), Some(y)) => x == y,
            _ => false,
        },
    }
}

/// Python's `max` over floats: the first of the largest.
pub fn py_max(xs: impl IntoIterator<Item = f64>) -> Option<f64> {
    xs.into_iter().fold(None, |m, x| match m {
        Some(b) if x.partial_cmp(&b) != Some(std::cmp::Ordering::Greater) => Some(b),
        _ => Some(x),
    })
}

/// Python's `min` over floats: the first of the smallest.
pub fn py_min(xs: impl IntoIterator<Item = f64>) -> Option<f64> {
    xs.into_iter().fold(None, |m, x| match m {
        Some(b) if x.partial_cmp(&b) != Some(std::cmp::Ordering::Less) => Some(b),
        _ => Some(x),
    })
}

/// The text of a file read as Python's `open(path).read()` reads it:
/// strict UTF-8, line endings turned into `\n`.
pub fn read_text(path: &Path) -> std::io::Result<String> {
    let bytes = std::fs::read(path)?;
    let text = String::from_utf8(bytes)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    Ok(if text.contains('\r') {
        text.replace("\r\n", "\n").replace('\r', "\n")
    } else {
        text
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paths() {
        assert_eq!(
            normpath(Path::new("/a/b/../c/./d")),
            PathBuf::from("/a/c/d")
        );
        assert_eq!(normpath(Path::new("/..")), PathBuf::from("/"));
        assert_eq!(normpath(Path::new("a/../..")), PathBuf::from(".."));
        assert_eq!(basename(Path::new("/x/trace_1.txt")), "trace_1.txt");
    }

    #[test]
    fn equality() {
        assert!(py_eq(&Value::Int(1), &Value::Float(1.0)));
        assert!(py_eq(&Value::Bool(true), &Value::Int(1)));
        assert!(!py_eq(&Value::Str("1".into()), &Value::Int(1)));
        assert!(!py_eq(&Value::Null, &Value::Int(0)));
        assert_eq!(py_max([0.0, -0.0]).map(f64::to_bits), Some(0f64.to_bits()));
        assert_eq!(
            py_max([-0.0, 0.0]).map(f64::to_bits),
            Some((-0f64).to_bits())
        );
    }
}
