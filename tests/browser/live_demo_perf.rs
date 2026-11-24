//! Map frame budget on a genome-scale synthetic document.

use super::*;
use std::io::Read;
use std::thread;

const MEGABASE: i128 = 1_000_000;
const CHROM: &str = "chr1";
const HEAD: &str = "function draw() {\n      raf = 0;";
const TIMING: &str = r#"(() => {
    const raf = window.requestAnimationFrame.bind(window);
    let pending = 0, start = 0;
    window.mapFrames = [];
    window.addEventListener('pointermove', () => { start = performance.now(); }, true);
    window.addEventListener('pointermove', () => {
        if (start) pending += performance.now() - start;
        start = 0;
    });
    window.requestAnimationFrame = callback => raf(t => {
        const begin = performance.now();
        callback(t);
        window.mapFrames.push(performance.now() - begin + pending);
        pending = 0;
    });
})()"#;

fn object(fields: impl IntoIterator<Item = (&'static str, Value)>) -> Value {
    Value::Object(fields.into_iter().map(|(k, v)| (k.into(), v)).collect())
}

fn fields(value: &mut Value) -> &mut indexmap::IndexMap<String, Value> {
    let Value::Object(fields) = value else {
        panic!("expected JSON object")
    };
    fields
}

fn text(value: &Value) -> &str {
    let Value::Str(value) = value else {
        panic!("expected text")
    };
    value
}

fn number(value: &Value) -> i128 {
    let Value::Int(value) = value else {
        panic!("expected integer")
    };
    *value
}

fn numbers(values: impl IntoIterator<Item = i128>) -> Value {
    Value::Array(values.into_iter().map(Value::Int).collect())
}

fn read_gzip(path: &Path) -> Value {
    let mut text = String::new();
    flate2::read::GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    loads(&text).unwrap()
}

fn write_gzip(path: &Path, value: &Value) {
    let mut out =
        flate2::write::GzEncoder::new(fs::File::create(path).unwrap(), flate2::Compression::fast());
    out.write_all(
        pycompat::json::dumps(value, &Default::default())
            .unwrap()
            .as_bytes(),
    )
    .unwrap();
    out.finish().unwrap();
}

fn fixture(map: &Path, genes: &Path) {
    let mut doc = read_gzip(map);
    let Some(Value::Array(chroms)) = prop(&doc, "ref").and_then(|o| prop(o, "chroms")) else {
        panic!("map has no chromosomes");
    };
    let lengths: Vec<_> = chroms
        .iter()
        .map(|c| (text(prop(c, "n").unwrap()), number(prop(c, "len").unwrap())))
        .collect();
    let chr1 = lengths[0].1;
    let total: i128 = lengths.iter().map(|(_, n)| n).sum();
    let mut run_count: i128 = 0;
    let mut runs = indexmap::IndexMap::new();
    let mut density = indexmap::IndexMap::from([("size".into(), Value::Int(MEGABASE))]);
    for (chrom, length) in lengths {
        let count = (4000 * length + total / 2) / total;
        let stride = length / (count + 1);
        let mut coverage = vec![0; (length / MEGABASE + 1) as usize];
        let mut starts = Vec::new();
        let mut spans = Vec::new();
        let mut prev = 0;
        for i in 1..=count {
            run_count += 1;
            let start = i * stride - stride / 2;
            let span = if run_count % 50 == 0 {
                6 * MEGABASE
            } else {
                100_000 + run_count * 7919 % 1_300_000
            };
            let end = (start + span).min(length - 1);
            starts.push(start - prev);
            spans.push(end - start);
            prev = start;
            for bin in start / MEGABASE..=(end - 1) / MEGABASE {
                coverage[bin as usize] += end.min((bin + 1) * MEGABASE) - start.max(bin * MEGABASE);
            }
        }
        runs.insert(
            chrom.to_owned(),
            object([
                ("ds", numbers(starts)),
                ("markers", numbers(spans.iter().map(|l| l / 1400 + 1))),
                (
                    "long",
                    numbers(spans.iter().map(|l| i128::from(*l >= 5 * MEGABASE))),
                ),
                ("l", numbers(spans)),
                ("q", numbers(std::iter::repeat_n(50, count as usize))),
                ("art", numbers(std::iter::repeat_n(0, count as usize))),
            ]),
        );
        density.insert(chrom.to_owned(), numbers(coverage));
    }
    let top = fields(&mut doc);
    let per = fields(fields(top.get_mut("per").unwrap()).get_mut(MAIN).unwrap());
    per.insert("roh".into(), Value::Object(runs));
    per.insert("roh_bins".into(), Value::Object(density));
    let roh = fields(
        fields(top.get_mut("layers").unwrap())
            .get_mut("roh")
            .unwrap(),
    );
    roh.insert("state".into(), Value::Str("ok".into()));
    roh.insert("count".into(), Value::Int(run_count));
    write_gzip(map, &doc);

    let mut deltas = Vec::new();
    let mut spans = Vec::new();
    let mut names = Vec::new();
    let mut types = Vec::new();
    let mut strands = Vec::new();
    let mut exons = Vec::new();
    let mut prev = 0;
    for i in 0..5000 {
        let start = 10_000 + i * (chr1 - 200_000) / 5000;
        let span = 20_000 + i * 104_729 % 100_000;
        deltas.push(start - prev);
        spans.push(span);
        prev = start;
        strands.push(Value::Str(if i % 2 == 0 { "-" } else { "+" }.into()));
        types.push(Value::Str(
            if i % 3 == 0 {
                "lncRNA"
            } else {
                "protein_coding"
            }
            .into(),
        ));
        names.push(Value::Str(format!("SYNG{}", i + 1)));
        for j in 0..8 {
            let exon = start + j * span / 8;
            exons.push((exon, exon + 300 + (i * 31 + j * 97) % 900));
        }
    }
    exons.sort_unstable();
    let mut merged: Vec<(i128, i128)> = Vec::new();
    for (start, end) in exons {
        if let Some(last) = merged.last_mut()
            && last.1 >= start
        {
            last.1 = last.1.max(end);
        } else {
            merged.push((start, end));
        }
    }
    let mut prev = 0;
    let mut exon_starts = Vec::new();
    let mut exon_lengths = Vec::new();
    for (start, end) in merged {
        exon_starts.push(start - prev);
        exon_lengths.push(end - start);
        prev = start;
    }
    write_gzip(
        genes,
        &object([
            ("v", Value::Int(1)),
            ("gencode", Value::Str("v50".into())),
            ("chrom", Value::Str(CHROM.into())),
            (
                "g",
                object([
                    ("ds", numbers(deltas)),
                    ("l", numbers(spans)),
                    ("strand", Value::Array(strands)),
                    ("name", Value::Array(names)),
                    ("type", Value::Array(types)),
                ]),
            ),
            (
                "x",
                object([("ds", numbers(exon_starts)), ("l", numbers(exon_lengths))]),
            ),
        ]),
    );
}

fn decimal(value: &Value) -> f64 {
    match value {
        Value::Int(n) => *n as f64,
        Value::Float(n) => *n,
        _ => panic!("expected browser number, got {value:?}"),
    }
}

fn mouse(b: &mut Browser, action: &str, x: f64, y: f64, down: u8) {
    call(
        b,
        "Input.dispatchMouseEvent",
        &format!(
            r#"{{"type":"{action}","x":{x},"y":{y},"button":"left","clickCount":1,"buttons":{down}}}"#
        ),
    );
}

fn pan(b: &mut Browser, base: &str) -> (Vec<f64>, f64, f64, f64) {
    go(b, &run_at(base, MAIN, ""));
    map_level(b, "genome");
    let parse = decimal(&eval(b, &format!("Number({MAP}.dataset.pgpMapParseMs)")));
    let paint = decimal(&eval(b, &format!("Number({MAP}.dataset.pgpMapPaintMs)")));
    let mut times = Vec::new();
    for (hash, level) in [
        ("#chr1:100.0-160.0M", "chromosome"),
        ("#chr1:100.0-105.0M", "region"),
        ("#chr1:100.000-100.200M", "gene"),
    ] {
        eval(b, &format!("location.hash = {}; true", js_str(hash)));
        map_level(b, level);
        if level != "chromosome" {
            require(b, &format!("Number({MAP}.dataset.genes) > 0"));
        }
        thread::sleep(Duration::from_millis(650));
        let Value::Array(point) = eval(
            b,
            r#"(() => {
            const canvas = document.querySelector('.gm-canvas');
            canvas.scrollIntoView({block: 'center', behavior: 'instant'});
            const r = canvas.getBoundingClientRect();
            return [r.left + r.width / 2, r.bottom - 30];
        })()"#,
        ) else {
            panic!("missing map canvas")
        };
        let (mut x, y) = (decimal(&point[0]), decimal(&point[1]));
        require(
            b,
            &format!(
                "document.elementFromPoint({x}, {y}) === document.querySelector('.gm-canvas')"
            ),
        );
        mouse(b, "mouseMoved", x, y, 0);
        eval(b, "window.mapFrames = []; true");
        mouse(b, "mousePressed", x, y, 1);
        for i in 0..100 {
            x += if (i / 25) % 2 == 0 { 4.0 } else { -4.0 };
            mouse(b, "mouseMoved", x, y, 1);
            thread::sleep(Duration::from_millis(20));
        }
        mouse(b, "mouseReleased", x, y, 0);
        map_level(b, level);
        let Value::Array(frames) = eval(b, "window.mapFrames") else {
            panic!("no animation frames")
        };
        times.extend(frames.iter().map(decimal));
    }
    let heap = decimal(&eval(b, "performance.memory.usedJSHeapSize")) / 1e6;
    times.sort_by(f64::total_cmp);
    (times, paint, parse, heap)
}

const SCROLL_TABLE: &str = r#"(async () => {
    const table = document.querySelector('.gm-inview-wrap');
    const keys = new Set();
    for (let i = 0; i < 400; i++) {
        for (const row of table.querySelectorAll('tbody tr[data-key]')) keys.add(row.dataset.key);
        if (table.scrollTop + table.clientHeight >= table.scrollHeight - 1) break;
        table.scrollTop += 240;
        await new Promise(done => requestAnimationFrame(() => requestAnimationFrame(done)));
    }
    return [...keys].sort();
})()"#;

fn table_window(b: &mut Browser, base: &str) -> bool {
    go(b, &run_at(base, MAIN, "#chr1"));
    map_level(b, "chromosome");
    require(
        b,
        "document.querySelector('.gm-canvas').pgpDrawn.length > 200",
    );
    css(b, ".gm-table tbody tr[data-key]", 12);
    require(
        b,
        "document.querySelector('.gm-inview-foot').textContent.includes('more in view')",
    );
    eval(b, "document.querySelector('.gm-showall').click(); true");
    require(
        b,
        r#"(() => {
        const wrap = document.querySelector('.gm-inview-wrap');
        const count = wrap.querySelectorAll('tbody tr[data-key]').length;
        return wrap.classList.contains('gm-virtual') && count > 12
            && count < document.querySelector('.gm-canvas').pgpDrawn.length;
    })()"#,
    );
    let seen = call(
        b,
        "Runtime.evaluate",
        &format!(
            r#"{{"expression":{},"awaitPromise":true,"returnByValue":true}}"#,
            js_str(SCROLL_TABLE)
        ),
    );
    let got = prop(&seen, "result").and_then(|v| prop(v, "value"));
    got == Some(&eval(
        b,
        "document.querySelector('.gm-canvas').pgpDrawn.slice().sort()",
    ))
}

pub(super) fn check(b: &mut Browser, base: &str, demo: &Path) {
    viewport(b, 1440, 900);
    let map = demo.join(format!("runs/{MAIN}/map.json"));
    let genes = demo.join("map/genes/chr1.json");
    let originals = [fs::read(&map).unwrap(), fs::read(&genes).unwrap()];
    fixture(&map, &genes);
    let response = call(
        b,
        "Page.addScriptToEvaluateOnNewDocument",
        &format!(r#"{{"source":{}}}"#, js_str(TIMING)),
    );
    let id = text(prop(&response, "identifier").unwrap());
    let (times, paint, parse, heap) = pan(b, base);
    assert!(
        times.len() >= 200,
        "pan drew {} frames, expected at least 200",
        times.len()
    );
    let median = times[times.len() / 2];
    let p95 = times[times.len() * 95 / 100];
    println!(
        "map pan: {} frames, median {median:.2} ms, p95 {p95:.2} ms, first paint {paint:.2} ms, parse {parse:.2} ms, heap {heap:.1} MB",
        times.len()
    );
    assert!(median <= 8.0 && p95 <= 16.0, "map frame budget exceeded");
    assert!(
        paint <= 150.0 && parse <= 30.0,
        "first map render exceeded its budget"
    );
    assert!(heap > 0.0 && heap <= 60.0, "map heap exceeded its budget");

    let script = demo.join("assets/map.js");
    let js = fs::read_to_string(&script).unwrap();
    assert!(
        table_window(b, base),
        "the virtual window missed a painted row"
    );
    let window = "var a = clamp(Math.floor(inWrap.scrollTop / ROW_PX) - 10, 0, cap)";
    assert_eq!(js.matches(window).count(), 1);
    fs::write(&script, js.replace(window, "var a = 0")).unwrap();
    assert!(
        !table_window(b, base),
        "the fixed first window control should fail"
    );
    assert_eq!(
        js.matches(HEAD).count(),
        1,
        "the slow draw control must patch exactly once"
    );
    let slow = js.replace(HEAD, "function draw() {\n      const until = performance.now() + 20;\n      while (performance.now() < until) { void 0; }\n      raf = 0;");
    fs::write(&script, slow).unwrap();
    let (times, _, _, _) = pan(b, base);
    assert!(!times.is_empty(), "the slow draw control drew nothing");
    assert!(
        times[times.len() * 95 / 100] > 16.0,
        "the slow draw control fit the frame budget"
    );
    fs::write(script, js).unwrap();
    for (path, bytes) in [map, genes].into_iter().zip(originals) {
        fs::write(path, bytes).unwrap();
    }
    call(
        b,
        "Page.removeScriptToEvaluateOnNewDocument",
        &format!(r#"{{"identifier":{}}}"#, js_str(id)),
    );
}
