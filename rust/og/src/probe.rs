//! `og pycompat-probe`: reads one case per line and writes one result per
//! line, for `tests/rust/pycompat-diff.sh` to compare with what CPython
//! writes for the same cases (`tests/rust/pycompat_oracle.py`).
//!
//! A case is tab-separated: the operation, then its arguments. Floats travel
//! as the 16 hex digits of their bits, text as the hex of its UTF-8, so no
//! value is rounded or mangled on the way. A result is the Python result,
//! or the name of the exception Python raised, escaped by [`escape`].

use og_io::{Lines, Utf8};
use pycompat::{fmt, json, num};
use std::io::{self, BufRead, Write};

/// `\` `\n` `\r` `|` and NUL escaped, so every result is one line of
/// text and a list of lines can be joined with `|`.
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '|' => out.push_str("\\|"),
            '\0' => out.push_str("\\0"),
            c => out.push(c),
        }
    }
    out
}

fn bad(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn float(s: &str) -> io::Result<f64> {
    u64::from_str_radix(s, 16)
        .map(f64::from_bits)
        .map_err(|_| bad(format!("bad float bits {s:?}")))
}

fn bytes(s: &str) -> io::Result<Vec<u8>> {
    if !s.len().is_multiple_of(2) {
        return Err(bad("odd hex length"));
    }
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).map_err(|_| bad("bad hex")))
        .collect()
}

fn text(s: &str) -> io::Result<String> {
    String::from_utf8(bytes(s)?).map_err(|_| bad("hex text is not UTF-8"))
}

fn int<T: std::str::FromStr>(s: &str) -> io::Result<T> {
    s.parse().map_err(|_| bad(format!("bad integer {s:?}")))
}

fn floats(s: &str) -> io::Result<Vec<f64>> {
    if s.is_empty() {
        Ok(Vec::new())
    } else {
        s.split(',').map(float).collect()
    }
}

fn dump_opts(s: &str) -> io::Result<json::Dump<'static>> {
    let f: Vec<&str> = s.split(',').collect();
    let [indent, sep, ascii, sort, nan] = f[..] else {
        return Err(bad(format!("bad json options {s:?}")));
    };
    Ok(json::Dump {
        indent: if indent == "-" {
            None
        } else {
            Some(int(indent)?)
        },
        separators: match sep {
            "d" => None,
            "c" => Some((",", ":")),
            "p" => Some((", ", ": ")),
            "i" => Some((",", ": ")),
            _ => return Err(bad(format!("bad separators {sep:?}"))),
        },
        ensure_ascii: ascii == "1",
        sort_keys: sort == "1",
        allow_nan: nan == "1",
    })
}

/// One case's result.
fn eval(f: &[&str]) -> io::Result<String> {
    let arg = |i: usize| {
        f.get(i)
            .copied()
            .ok_or_else(|| bad(format!("{}: missing argument {i}", f[0])))
    };
    Ok(match f[0] {
        "repr" => fmt::repr(float(arg(1)?)?),
        "f" => fmt::fixed(float(arg(1)?)?, int(arg(2)?)?),
        "e" => fmt::sci(float(arg(1)?)?, int(arg(2)?)?),
        "g" => fmt::general(float(arg(1)?)?, int(arg(2)?)?),
        "pct" => fmt::percent(float(arg(1)?)?, int(arg(2)?)?),
        "grpf" => fmt::fixed_grouped(float(arg(1)?)?, int(arg(2)?)?),
        "grp" => fmt::int_grouped(int(arg(1)?)?),
        "zpad" => fmt::zero_pad(int(arg(1)?)?, int(arg(2)?)?),
        "ljust" => fmt::ljust(&text(arg(1)?)?, int(arg(2)?)?),
        "rjust" => fmt::rjust(&text(arg(1)?)?, int(arg(2)?)?),
        "round0" => match fmt::round_int(float(arg(1)?)?) {
            Ok(v) => fmt::int_str(v),
            Err(e) => e.exception().to_owned(),
        },
        "round" => match fmt::round(float(arg(1)?)?, int(arg(2)?)?) {
            Ok(v) => fmt::repr(v),
            Err(e) => e.exception().to_owned(),
        },
        "float" => match num::parse_float(&text(arg(1)?)?) {
            Ok(v) => fmt::repr(v),
            Err(e) => e.exception().to_owned(),
        },
        "int" => match num::parse_int(&text(arg(1)?)?) {
            Ok(v) => v.to_string(),
            Err(e) => e.exception().to_owned(),
        },
        "fsum" => match num::fsum(floats(arg(1)?)?) {
            Ok(v) => fmt::repr(v),
            Err(e) => e.exception().to_owned(),
        },
        "sum" => fmt::repr(num::sum(floats(arg(1)?)?)),
        "json" => {
            let opts = dump_opts(arg(2)?)?;
            match json::loads(&text(arg(1)?)?).and_then(|v| json::dumps(&v, &opts)) {
                Ok(s) => s,
                Err(e) => e.exception().to_owned(),
            }
        }
        "lines" => {
            let lines: Vec<String> = Lines::new(io::Cursor::new(bytes(arg(1)?)?), Utf8::Replace)
                .collect::<io::Result<_>>()?;
            return Ok(lines
                .iter()
                .map(|l| escape(l))
                .collect::<Vec<_>>()
                .join("|"));
        }
        op => return Err(bad(format!("unknown operation {op:?}"))),
    })
}

pub fn run(input: impl BufRead, output: impl Write) -> io::Result<()> {
    let mut out = io::BufWriter::new(output);
    for line in input.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let result = eval(&fields)?;
        // A lines result is escaped line by line already.
        if fields[0] == "lines" {
            writeln!(out, "{result}")?;
        } else {
            writeln!(out, "{}", escape(&result))?;
        }
    }
    out.flush()
}
