//! Strict checks over raw Nextflow traces. Only one explicit trace path is accepted.
use pycompat::csv::{Quoting, Reader};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

const GB: u64 = 1 << 30;
const H: u64 = 3_600_000;
type Record = BTreeMap<String, String>;

fn records(path: &str, mut f: impl FnMut(Record) -> Result<(), String>) -> Result<(), String> {
    if path.contains(['*', '?', '[', ']']) {
        return Err("trace glob refused: pass exactly one explicit trace path".into());
    }
    let lines = super::smoke::csv_lines(path)?;
    let mut csv = Reader::new(lines, '\t', Quoting::Minimal);
    let header = csv
        .next()
        .ok_or("no trace header")?
        .map_err(|e| e.to_string())?;
    for row in csv {
        let row = row.map_err(|e| e.to_string())?;
        if row.is_empty() {
            continue;
        }
        if row.len() != header.len() {
            return Err(format!("{path}: row width does not match header"));
        }
        f(header.iter().cloned().zip(row).collect())?;
    }
    Ok(())
}
fn cell<'a>(r: &'a Record, k: &str) -> Result<&'a str, String> {
    r.get(k)
        .map(String::as_str)
        .ok_or_else(|| format!("missing column {k}"))
}
fn num(r: &Record, k: &str) -> Result<u64, String> {
    cell(r, k)?
        .parse()
        .map_err(|_| format!("column {k} is not an unsigned integer"))
}
fn traces(path: &str, out: &mut dyn Write) -> Result<Vec<Record>, String> {
    let mut rows: Vec<Record> = Vec::new();
    let mut ids = BTreeMap::new();
    records(path, |r| {
        let id = cell(&r, "task_id")?.to_owned();
        if let Some(&i) = ids.get(&id) {
            let first: &Record = &rows[i];
            let a = cell(first, "exit")?;
            let b = cell(&r, "exit")?;
            if (a != "-" && b != "-" && a != b)
                || r.iter()
                    .any(|(k, v)| k != "exit" && first.get(k) != Some(v))
            {
                return Err(format!("task {id} has two different rows in the trace"));
            }
            writeln!(
                out,
                "{}: attempt {} is written twice in the trace; counted once",
                cell(&r, "name")?,
                cell(&r, "attempt")?
            )
            .map_err(|e| e.to_string())?;
        } else {
            ids.insert(id, rows.len());
            rows.push(r);
        }
        Ok(())
    })?;
    if rows.is_empty() {
        return Err("no trace rows".into());
    }
    Ok(rows)
}
fn options(args: &[String]) -> Result<BTreeMap<&str, &str>, String> {
    if !args.len().is_multiple_of(2) {
        return Err("expected named option and value; pass exactly one trace path".into());
    }
    let mut opts = BTreeMap::new();
    for pair in args.as_chunks::<2>().0 {
        if !pair[0].starts_with("--") || opts.insert(pair[0].as_str(), pair[1].as_str()).is_some() {
            return Err("unexpected or repeated argument".into());
        }
    }
    Ok(opts)
}
fn opt<'a>(opts: &BTreeMap<&str, &'a str>, key: &str) -> Result<&'a str, String> {
    opts.get(key)
        .copied()
        .ok_or_else(|| format!("{key} is required"))
}
fn glob(pattern: &str, text: &str) -> bool {
    // Escaping makes punctuation literal. Only '*' and '?' have glob meaning.
    let re = regex::escape(pattern)
        .replace("\\*", ".*")
        .replace("\\?", ".");
    regex::Regex::new(&format!("(?s)^{re}$")).is_ok_and(|r| r.is_match(text))
}
fn attempts(path: &str, opts: &BTreeMap<&str, &str>, out: &mut dyn Write) -> Result<bool, String> {
    for k in opts.keys() {
        if !["--name-glob", "--want", "--nextflow-exit"].contains(k) {
            return Err(format!("unknown option {k}"));
        }
    }
    let pattern = opt(opts, "--name-glob")?;
    let want: Vec<(u64, String, u64)> = opt(opts, "--want")?
        .split(',')
        .map(|t| {
            let parts: Vec<_> = t.split(':').collect();
            if parts.len() != 3 {
                return Err("--want must be attempt:exit:GB,...".into());
            }
            Ok((
                parts[0].parse().map_err(|_| "bad attempt")?,
                parts[1].to_owned(),
                parts[2].parse().map_err(|_| "bad memory")?,
            ))
        })
        .collect::<Result<_, String>>()?;
    let mut got = Vec::new();
    for r in traces(path, out)? {
        if glob(pattern, cell(&r, "name")?) {
            got.push((
                num(&r, "attempt")?,
                cell(&r, "exit")?.to_owned(),
                num(&r, "memory")? / GB,
            ));
        }
    }
    got.sort();
    let tuples = |xs: &[(u64, String, u64)]| {
        format!(
            "[{}]",
            xs.iter()
                .map(|(a, e, m)| format!("({a}, '{e}', {m})"))
                .collect::<Vec<_>>()
                .join(", ")
        )
    };
    let label = if pattern.contains("ALIGN_MINIMAP2") {
        "ALIGN_MINIMAP2"
    } else {
        pattern
    };
    writeln!(
        out,
        "{label} attempts (attempt, exit, GB): {}; want {}",
        tuples(&got),
        tuples(&want)
    )
    .map_err(|e| e.to_string())?;
    let ok = got == want && opt(opts, "--nextflow-exit")? == "0";
    writeln!(
        out,
        "{}",
        if ok {
            "OK"
        } else {
            "FAIL: a memory kill did not reach Nextflow as exit 137 and a retry that passes"
        }
    )
    .map_err(|e| e.to_string())?;
    Ok(ok)
}
fn requests(path: &str, opts: &BTreeMap<&str, &str>, out: &mut dyn Write) -> Result<bool, String> {
    for k in opts.keys() {
        if !["--want", "--executor-cpus", "--executor-mem-gb"].contains(k) {
            return Err(format!("unknown option {k}"));
        }
    }
    let mut want = BTreeMap::new();
    records(opt(opts, "--want")?, |r| {
        let name = cell(&r, "process")?.to_owned();
        let v = (
            num(&r, "cpus")?,
            num(&r, "memory_gb")?
                .checked_mul(GB)
                .ok_or("memory overflow")?,
            num(&r, "time_h")?.checked_mul(H).ok_or("time overflow")?,
        );
        if want.insert(name, v).is_some() {
            return Err("duplicate process in --want".into());
        }
        Ok(())
    })?;
    if want.is_empty() {
        return Err("empty request expectations".into());
    }
    let rows = traces(path, out)?;
    type Seen<'a> = BTreeMap<&'a str, BTreeMap<u64, BTreeMap<&'a str, Vec<(u64, u64, u64)>>>>;
    let mut seen: Seen<'_> = BTreeMap::new();
    let mut all_cpus = 0u64;
    let mut all_mem = 0u64;
    let mut bad = 0;
    for r in &rows {
        let cpus = num(r, "cpus")?;
        let mem = num(r, "memory")?;
        all_cpus = all_cpus.checked_add(cpus).ok_or("CPU sum overflow")?;
        all_mem = all_mem.checked_add(mem).ok_or("memory sum overflow")?;
        let name = cell(r, "name")?;
        let proc = name
            .split(" (")
            .next()
            .unwrap_or(name)
            .rsplit(':')
            .next()
            .unwrap_or(name);
        if !want.contains_key(proc) {
            continue;
        }
        let attempt = num(r, "attempt")?;
        let time = num(r, "time")?;
        writeln!(
            out,
            "{name}: attempt {attempt}, exit {}, cpus {cpus}, memory {} GB, time {} h",
            cell(r, "exit")?,
            mem / GB,
            time / H
        )
        .map_err(|e| e.to_string())?;
        seen.entry(proc)
            .or_default()
            .entry(attempt)
            .or_default()
            .entry(name)
            .or_default()
            .push((cpus, mem, time));
    }
    for (proc, &expected) in &want {
        let attempts = seen.get(proc.as_str());
        if attempts.is_some_and(|a| a.keys().any(|k| ![1, 2].contains(k))) {
            writeln!(out, "FAIL {proc}: want only a first attempt and one retry")
                .map_err(|e| e.to_string())?;
            bad += 1;
        }
        let first = attempts.and_then(|a| a.get(&1));
        let second = attempts.and_then(|a| a.get(&2));
        if first.is_none_or(|a| a.is_empty()) {
            writeln!(out, "FAIL {proc}: no first attempt in the trace")
                .map_err(|e| e.to_string())?;
            bad += 1;
            continue;
        }
        let names: BTreeSet<_> = first
            .into_iter()
            .chain(second)
            .flat_map(|a| a.keys())
            .collect();
        for name in names {
            let a = first.and_then(|a| a.get(name)).map_or(0, Vec::len);
            let b = second.and_then(|a| a.get(name)).map_or(0, Vec::len);
            if a != 1 || b != 1 {
                writeln!(
                    out,
                    "FAIL {name}: {a} first and {b} second attempts in the trace, want one of each"
                )
                .map_err(|e| e.to_string())?;
                bad += 1;
            }
        }
        for got in first.into_iter().flat_map(|a| a.values()).flatten() {
            if *got != expected {
                writeln!(out,"FAIL {proc}: first attempt asks for cpus {}, memory {} GB, time {} h; want {}, {} GB, {} h",got.0,got.1/GB,got.2/H,expected.0,expected.1/GB,expected.2/H).map_err(|e| e.to_string())?;
                bad += 1;
            }
        }
        for got in second.into_iter().flat_map(|a| a.values()).flatten() {
            if Some(got.1) != expected.1.checked_mul(2) {
                writeln!(
                    out,
                    "FAIL {proc}: the retry asks for {} GB, want {} GB",
                    got.1 / GB,
                    expected.1 / GB * 2
                )
                .map_err(|e| e.to_string())?;
                bad += 1;
            }
        }
    }
    let cpus: u64 = opt(opts, "--executor-cpus")?
        .parse()
        .map_err(|_| "invalid executor CPUs")?;
    let mem: u64 = opt(opts, "--executor-mem-gb")?
        .parse()
        .map_err(|_| "invalid executor memory")?;
    writeln!(out,"{} task attempts ask for {all_cpus} CPUs and {} GB added together; the executor is told it has {cpus} CPUs and {mem} GB",rows.len(),all_mem/GB).map_err(|e| e.to_string())?;
    if all_cpus > cpus || Some(all_mem) > mem.checked_mul(GB) {
        writeln!(out,"FAIL: executor.cpus or executor.memory of this step is below that sum; raise exec_cpus / exec_mem_gb").map_err(|e| e.to_string())?;
        bad += 1;
    }
    writeln!(out, "{bad} failure(s)").map_err(|e| e.to_string())?;
    Ok(bad == 0)
}

/// Checks `requests TRACE ...` or `attempts TRACE ...`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let result = (|| {
        if args.len() < 2 {
            return Err("usage: og dev nf-trace requests|attempts TRACE [options]".into());
        }
        let opts = options(&args[2..])?;
        match args[0].as_str() {
            "requests" => requests(&args[1], &opts, out),
            "attempts" => attempts(&args[1], &opts, out),
            _ => Err("unknown nf-trace verb".into()),
        }
    })();
    match result {
        Ok(true) => 0,
        Ok(false) => 1,
        Err(e) => {
            let _ = writeln!(err, "FAIL: {e}");
            2
        }
    }
}
