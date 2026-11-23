//! Build the public live site from the invented sample and the live server's
//! own state, results, family pages, and policy.

use crate::family::{self, inside, realpath};
use crate::map::{Ended, OG_SECONDS, run_og};
use crate::runs::{self, Run};
use crate::server::{AccessLog, Config};
use crate::state::{Live, dumps_compact, get, obj, s};
use crate::static_files as pages;
use crate::{csp, ingest, plugins, results, state};
use jiff::Timestamp;
use og_dev::demo::{RUN_STARTED, Scratch, find_repo, put, read, utc};
use og_dev::versions_env;
use og_runlog::host::{FsStat, Host};
use og_runlog::steps::read_steps;
use pycompat::json::{Value, loads};
use std::collections::BTreeSet;
use std::ffi::OsStr;
use std::fs;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const SAMPLE: &str = "DEMO-001";
const FAMILY: &str = "DEMO-FAMILY";
const OTHERS: [(&str, &str, &str); 2] = [("DEMO-002", "female", "2"), ("DEMO-003", "male", "3")];
const FIXTURES: &str = "tests/fixtures/live";
const USAGE: &str = "usage: og live demo --out DIR [--repo DIR] [--genome-dir DIR] [--records N] [--work DIR] [--og OG]";
const INTERVAL_MS: i128 = 600;
const GIB: u64 = 1 << 30;

struct Args {
    output: PathBuf,
    repo: PathBuf,
    genome: Option<PathBuf>,
    records: usize,
    work: Option<PathBuf>,
    og: PathBuf,
}

impl Args {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut output = None;
        let mut repo = None;
        let mut genome = None;
        let mut records = 20_000;
        let mut work = None;
        let mut og = None;
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            let value = args
                .next()
                .filter(|v| !v.starts_with("--"))
                .ok_or_else(|| format!("{flag} takes a value"))?;
            match flag.as_str() {
                "--out" => output = Some(PathBuf::from(value)),
                "--repo" => repo = Some(PathBuf::from(value)),
                "--genome-dir" => genome = Some(realpath(Path::new(value))),
                "--records" => {
                    records = value
                        .parse()
                        .map_err(|_| "--records takes a non-negative integer")?
                }
                "--work" => work = Some(realpath(Path::new(value))),
                "--og" => og = Some(PathBuf::from(value)),
                _ => return Err(format!("unknown argument {flag}")),
            }
        }
        let output = realpath(&output.ok_or("--out is required")?);
        let repo = repo
            .or_else(find_repo)
            .or_else(|| {
                let exe = std::env::current_exe().ok()?;
                exe.ancestors()
                    .find(|p| p.join("assets/steps.tsv").is_file())
                    .map(Path::to_path_buf)
            })
            .ok_or("cannot find the checkout; give --repo")?;
        let repo = realpath(&repo);
        if repo.starts_with(&output)
            || work.as_ref().is_some_and(|p| overlaps(p, &output))
            || genome.as_ref().is_some_and(|p| {
                overlaps(p, &output) || work.as_ref().is_some_and(|w| overlaps(p, w))
            })
        {
            return Err("--out, --work, the checkout and the input genome directory must not overwrite one another".into());
        }
        let og = og
            .or_else(|| std::env::var_os("OG").map(PathBuf::from))
            .or_else(|| std::env::current_exe().ok())
            .ok_or("cannot find og; give --og")?;
        Ok(Self {
            output,
            repo,
            genome,
            records,
            work,
            og,
        })
    }
}

fn overlaps(a: &Path, b: &Path) -> bool {
    a.starts_with(b) || b.starts_with(a)
}

fn path_text(path: &Path) -> Result<&str, String> {
    path.to_str()
        .ok_or_else(|| format!("{} is not a UTF-8 path", path.display()))
}

fn io_at<T>(path: &Path, result: io::Result<T>) -> Result<T, String> {
    result.map_err(|e| format!("{}: {e}", path.display()))
}

fn copy_tree(src: &Path, dst: &Path) -> Result<(), String> {
    fn copy_inside(src: &Path, dst: &Path, root: &Path) -> Result<(), String> {
        io_at(dst, fs::create_dir_all(dst))?;
        for entry in io_at(src, fs::read_dir(src))? {
            let entry = io_at(src, entry)?;
            let from = entry.path();
            let to = dst.join(entry.file_name());
            // An optional input is an invented sample folder. Never import a
            // symlink into another run or follow a symlinked directory.
            let kind = io_at(&from, entry.file_type())?;
            if kind.is_dir() {
                copy_inside(&from, &to, root)?;
            } else if (kind.is_file() || kind.is_symlink()) && inside(&from, root) && from.is_file()
            {
                io_at(&to, fs::copy(&from, &to))?;
            }
        }
        Ok(())
    }
    copy_inside(src, dst, &realpath(src))
}

fn og(args: &Args, words: &[&str], input: &[u8]) -> Result<Vec<u8>, String> {
    let words: Vec<_> = words.iter().map(|w| (*w).to_owned()).collect();
    let (ended, bytes, stderr) = run_og(
        &args.og,
        &words,
        input,
        std::time::Duration::from_secs(OG_SECONDS),
    );
    if ended == Ended::Code(0) {
        Ok(bytes)
    } else {
        Err(format!("og {}: {ended:?}: {stderr}", words[..2].join(" ")))
    }
}

fn json_file(path: &Path) -> Result<Value, String> {
    loads(&read(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

fn put_json(path: &Path, value: &Value) -> Result<(), String> {
    put(path, &dumps_compact(value))
}

fn items(value: Option<&Value>) -> &[Value] {
    match value {
        Some(Value::Array(items)) => items,
        _ => &[],
    }
}

fn string(value: Option<&Value>) -> Option<&str> {
    match value {
        Some(Value::Str(value)) => Some(value),
        _ => None,
    }
}

fn fields(value: &mut Value) -> Result<&mut state::Obj, String> {
    match value {
        Value::Object(fields) => Ok(fields),
        _ => Err("expected a JSON object".into()),
    }
}

/// One machine across every build, independent of the available disk and RAM.
struct DemoHost;

impl Host for DemoHost {
    fn statvfs(&self, _: &Path) -> io::Result<FsStat> {
        Ok(FsStat {
            frsize: 4096,
            bsize: 4096,
            blocks: 512 << 20,
            bfree: 230 << 20,
            bavail: 230 << 20,
            files: 128 << 20,
            ffree: 100 << 20,
        })
    }

    fn dev(&self, _: &Path) -> io::Result<u64> {
        Ok(1)
    }

    fn exists(&self, path: &Path) -> bool {
        path.exists()
    }

    fn physical_memory(&self) -> Option<u128> {
        Some(64 * u128::from(GIB))
    }
}

struct Demo<'a> {
    args: &'a Args,
    live: Live,
    start: Timestamp,
    end: Timestamp,
}

impl Demo<'_> {
    fn path(&self, rel: &str) -> PathBuf {
        self.args.repo.join(FIXTURES).join(rel)
    }

    fn sample(&self, name: &str, sex: &str, seed: &str, records: usize) -> Result<(), String> {
        let folder = self.live.root.join(name);
        let source = self.args.genome.as_ref().map(|r| r.join(name));
        if let Some(src) = source.filter(|p| p.is_dir()) {
            copy_tree(&src, &folder)?;
        } else {
            if name == SAMPLE && self.args.genome.is_some() {
                return Err("--genome-dir does not contain DEMO-001".into());
            }
            og(
                self.args,
                &[
                    "dev",
                    "demo-sample",
                    "--genome-dir",
                    path_text(&self.live.root)?,
                    "--repo",
                    path_text(&self.args.repo)?,
                    "--sample",
                    name,
                    "--sex",
                    sex,
                    "--seed",
                    seed,
                    "--records",
                    &records.to_string(),
                ],
                b"",
            )?;
        }

        let status_path = folder.join("logs/run_status.tsv");
        // Keep the sample's completed rows while collecting the finished
        // report. The replay gets its own plan with those results held back.
        let status = og_runlog::files::read_run_status(&status_path)
            .map_err(|e| format!("{}: {e}", status_path.display()))?
            .unwrap_or_default();
        let started = self.start.as_second();
        let ended = if name == SAMPLE {
            self.end.as_second()
        } else {
            started + 30 * 3600
        };
        let mut status_text = format!(
            "meta\tstarted_epoch\t{started}\nmeta\tstarted_utc\t{}\nmeta\tfinished_utc\t{}\nmeta\tdeclared_sex\t{sex}\n",
            utc(started),
            utc(ended)
        );
        for (step, outcome) in status.steps {
            status_text.push_str(&format!("step\t{step}\t{outcome}\n"));
        }
        put(&status_path, &status_text)?;
        self.manifest(&folder, sex, ended)?;
        og(
            self.args,
            &[
                "report",
                "render",
                "--sample",
                name,
                "--sample-dir",
                path_text(&folder)?,
                "--assets",
                path_text(&self.live.cfg.repo_assets)?,
                "--json",
                path_text(&folder.join("summary.json"))?,
                "-o",
                path_text(&folder.join(format!("{name}_report.html")))?,
                "--charts-json",
                path_text(&folder.join(format!("{name}_charts.json")))?,
                "--search-index",
                path_text(&folder.join(format!("{name}_search_index.json")))?,
            ],
            b"",
        )?;
        Ok(())
    }

    fn manifest(&self, folder: &Path, sex: &str, ended: i64) -> Result<(), String> {
        let root = self.args.genome.as_deref().unwrap_or(&self.live.root);
        let clinvar = root.join("clinvar/clinvar.vcf.gz");
        let date = og_io::open(&clinvar).ok().and_then(|reader| {
            reader
                .lines()
                .map_while(Result::ok)
                .take_while(|line| line.starts_with('#'))
                .find_map(|line| line.strip_prefix("##fileDate=").map(str::to_owned))
        });
        let hla = read(&self.live.root.join("t1k_idx/hlaidx/hla.dat"))
            .or_else(|_| read(&root.join("t1k_idx/hlaidx/hla.dat")))
            .ok()
            .and_then(|text| {
                text.lines().find_map(|line| {
                    line.split_once("(Rel. ")?
                        .1
                        .split_once(',')
                        .map(|(version, _)| format!("IPD-IMGT/HLA {version}"))
                })
            });
        let pins = versions_env::read(&self.args.repo.join("versions.env"))?;
        let bundle = versions_env::get(&pins, "PCGR_DATA_BUNDLE")
            .ok_or("versions.env has no PCGR_DATA_BUNDLE")?;
        let commit = Command::new("git")
            .args(["rev-parse", "HEAD"])
            .current_dir(&self.args.repo)
            .output()
            .ok()
            .filter(|r| r.status.success())
            .and_then(|r| String::from_utf8(r.stdout).ok());
        let mut manifest = format!(
            "run\twritten_utc\t{}\nrun\tstarted_utc\t{}\nrun\twritten_by\tog live demo\n\
             run\tdemo\tinvented data\nrun\tdeclared_sex\t{sex}\nrun\tgit_commit\t{}\n\
             data\tclinvar_file_date\t{}\ndata\tcpsr_clinvar_release\t{}\ndata\thla_database\t{}\n",
            utc(ended),
            utc(self.start.as_second()),
            commit.as_deref().map_or("source archive (demo)", str::trim),
            date.as_deref().unwrap_or("invented demo file"),
            format_args!("simulated (PCGR bundle {})", bundle.value),
            hla.as_deref().unwrap_or("simulated database")
        );
        for pin in pins {
            manifest.push_str(&format!("version\t{}\t{}\n", pin.name, pin.value));
        }
        put(&folder.join("run_manifest.tsv"), &manifest)
    }

    fn genes(&self) -> Result<(), String> {
        let src = self
            .args
            .repo
            .join("tests/fixtures/map/gencode-slice/reference");
        let dst = self.live.root.join("reference");
        copy_tree(&src, &dst)?;
        let gtf = io_at(&dst, fs::read_dir(&dst))?
            .filter_map(Result::ok)
            .map(|e| e.path())
            .find(|p| p.extension() == Some(OsStr::new("gtf")))
            .ok_or("the demo GENCODE slice has no GTF")?;
        let mut text = read(&gtf)?;
        // The synthetic ClinVar rows supply their own invented gene names
        // and spans. Add only those names, leaving the real GENCODE slice intact.
        let hits = read(
            &self
                .live
                .root
                .join(SAMPLE)
                .join(format!("clinvar/{SAMPLE}_clinvar_hits.vcf")),
        )?;
        let mut seen = BTreeSet::new();
        for line in hits.lines().filter(|l| !l.starts_with('#')) {
            let cells: Vec<_> = line.split('\t').collect();
            let [chrom, pos, _, reference, _, _, _, info, ..] = cells.as_slice() else {
                continue;
            };
            let Some(gene) = info
                .split(';')
                .find_map(|v| v.strip_prefix("GENEINFO="))
                .and_then(|v| v.split_once(':'))
                .map(|(g, _)| g)
                .filter(|g| matches!(*g, "DEMOA" | "DEMOB" | "DEMOC"))
            else {
                continue;
            };
            if seen.insert(gene) {
                let begin = pos
                    .parse::<u64>()
                    .map_err(|_| "invalid invented finding position")?;
                let end = begin + reference.len().saturating_sub(1) as u64;
                text.push_str(&format!(
                    "{chrom}\tOwnGenome\tgene\t{begin}\t{end}\t.\t+\t.\tgene_id \"{gene}\"; gene_name \"{gene}\"; gene_type \"protein_coding\";\n"
                ));
            }
        }
        put(&gtf, &text)
    }

    fn recheck(&self, checked: Timestamp) -> Result<(), String> {
        let mut check = json_file(&self.path("demo/DEMO-001/recheck/state.json"))?;
        let check_o = fields(&mut check)?;
        check_o.insert("checked_utc".into(), s(&checked.to_string()));
        let manifest = og_runlog::files::read_manifest(Some(
            &self.live.root.join(SAMPLE).join("run_manifest.tsv"),
        ))
        .map_err(|e| e.to_string())?;
        if let Some(Value::Array(items)) = check_o.get_mut("items") {
            for item in items {
                let item = fields(item)?;
                let key = string(item.get("item"))
                    .ok_or("re-check item has no key")?
                    .to_owned();
                let affected = if string(item.get("class")) == Some("no-newer-release") {
                    Vec::new()
                } else {
                    og(
                        self.args,
                        &[
                            "recheck",
                            "--affected",
                            &key,
                            "--steps",
                            path_text(&self.live.cfg.steps)?,
                        ],
                        b"",
                    )?
                };
                let source = if key == "clinvar_file_date" {
                    "data"
                } else {
                    "version"
                };
                if let Some(pin) = manifest.get(&(source.to_owned(), key.clone())) {
                    item.insert("ran_with".into(), s(pin));
                    item.insert("installed".into(), s(pin));
                    if source == "version" {
                        item.insert("upstream".into(), s(pin));
                    }
                }
                item.insert(
                    "steps".into(),
                    Value::Array(
                        String::from_utf8_lossy(&affected)
                            .split_whitespace()
                            .map(s)
                            .collect(),
                    ),
                );
            }
        }
        put_json(
            &self.live.root.join(SAMPLE).join("recheck/state.json"),
            &check,
        )
    }

    fn replay(&self, events: &[(Timestamp, Value)]) -> Result<Vec<Value>, String> {
        let run_folder = self.live.root.join(SAMPLE);
        let mut held = Vec::new();
        for file in [
            format!("{SAMPLE}_report.html"),
            format!("{SAMPLE}_charts.json"),
            format!("{SAMPLE}_search_index.json"),
            "summary.json".into(),
        ] {
            let src = run_folder.join(file);
            let dst = src.with_extension("held");
            io_at(&src, fs::rename(&src, &dst))?;
            held.push((src, dst));
        }
        let check = run_folder.join("recheck/state.json");
        if check.is_file() {
            io_at(&check, fs::remove_file(&check))?;
        }
        let selected: BTreeSet<_> = events
            .iter()
            .filter_map(|(_, event)| state::task_name(string(get(get(event, "trace")?, "name"))?).0)
            .collect();
        let mut plan = format!(
            "meta\tstarted_epoch\t{}\nmeta\tstarted_utc\t{}\nmeta\tdeclared_sex\tmale\n",
            self.start.as_second(),
            utc(self.start.as_second())
        );
        for step in &self.live.steps {
            if step.processes.iter().any(|p| selected.contains(p)) {
                plan.push_str(&format!("step\t{}\tplanned\n", step.step));
            }
        }
        let status = run_folder.join("logs/run_status.tsv");
        put(&status, &plan)?;
        let run = self
            .live
            .run(SAMPLE)
            .ok_or("the demo run was not discovered")?;
        let mut frames = Vec::new();
        let mut last = None;
        for (i, (time, event)) in events.iter().enumerate() {
            ingest::append(&self.live, SAMPLE, event).map_err(|e| e.to_string())?;
            if string(get(event, "event")) == Some("process_completed")
                && get(event, "trace")
                    .and_then(|t| string(get(t, "name")))
                    .is_some_and(|name| state::task_name(name).0.as_deref() == Some("HTML_REPORT"))
            {
                for (src, dst) in held.drain(..) {
                    io_at(&src, fs::rename(&dst, &src))?;
                }
            }
            if events.get(i + 1).is_some_and(|(next, _)| next == time) {
                continue;
            }
            if i + 1 == events.len() {
                plan.push_str(&format!(
                    "meta\tfinished_utc\t{}\n",
                    utc(self.end.as_second())
                ));
                put(&status, &plan)?;
            }
            let mut snapshot = self
                .live
                .replay_state(&run, time.as_second() as f64, &DemoHost);
            if frames.is_empty() && fault("host-path") {
                fields(&mut snapshot)?
                    .insert("planted_path".into(), s(path_text(&self.live.root)?));
            }
            let key = frame_key(&snapshot);
            if last.as_ref() != Some(&key) {
                last = Some(key);
                frames.push(frame(*time, snapshot));
            }
        }
        if !held.is_empty() {
            return Err("the demo weblog never completes HTML_REPORT".into());
        }
        let checked = self
            .end
            .checked_add(jiff::Span::new().hours(1))
            .map_err(|e| e.to_string())?;
        self.recheck(checked)?;
        frames.push(frame(
            checked,
            self.live
                .replay_state(&run, checked.as_second() as f64, &DemoHost),
        ));
        Ok(frames)
    }
}

fn frame_key(snapshot: &Value) -> Value {
    let done = items(get(snapshot, "steps"))
        .iter()
        .filter_map(|st| {
            matches!(string(get(st, "state")), Some("done" | "failed"))
                .then_some(get(st, "step"))
                .flatten()
                .cloned()
        })
        .collect();
    let mut plugins = get(snapshot, "plugins").cloned().unwrap_or(Value::Null);
    if let Value::Object(o) = &mut plugins {
        o.shift_remove("ops");
    }
    Value::Array(vec![
        get(snapshot, "run")
            .and_then(|r| get(r, "state"))
            .cloned()
            .unwrap_or(Value::Null),
        Value::Array(done),
        plugins,
    ])
}

fn frame(time: Timestamp, mut state: Value) -> Value {
    if let Value::Object(fields) = &mut state {
        let time = time.to_string();
        fields.insert("generated_utc".into(), s(&time));
        if let Some(plugins) = fields.get_mut("plugins") {
            relink(plugins, "../../");
        }
    }
    obj([("t", s(&time.to_string())), ("state", state)])
}

fn relink(value: &mut Value, prefix: &str) {
    match value {
        Value::Object(o) => o.values_mut().for_each(|v| relink(v, prefix)),
        Value::Array(a) => a.iter_mut().for_each(|v| relink(v, prefix)),
        Value::Str(v) if v.starts_with("/runs/") || v.starts_with("/family/") => {
            *v = format!("{prefix}{}", &v[1..]);
        }
        _ => {}
    }
}

fn meta(policy: &str) -> String {
    let policy: Vec<_> = policy
        .split("; ")
        .filter(|part| !part.starts_with("frame-ancestors") && !part.starts_with("sandbox"))
        .collect();
    format!(
        r#"<meta http-equiv="Content-Security-Policy" content="{}">"#,
        pages::esc(&policy.join("; "))
    )
}

fn page(html: String, prefix: &str) -> String {
    let theme = format!("{prefix}assets");
    html.replace("=\"/assets", &format!("=\"{theme}"))
        .replace("href=\"/\"", &format!("href=\"{prefix}\""))
        .replace("data-mode=\"live\"", "data-mode=\"replay\"")
        .replacen(
            "<meta charset=\"utf-8\">",
            &format!("<meta charset=\"utf-8\">\n{}", meta(csp::CSP)),
            1,
        )
}

fn copy_run_file(src: &Path, dst: &Path, name: &str) -> Result<(), String> {
    let warning = format!("Invented data: {name} is not a real person; nothing here is a result");
    match dst.extension().and_then(OsStr::to_str) {
        Some("txt") => put(dst, &format!("{warning}\n{}", read(src)?)),
        Some("html" | "htm") => {
            let mut html = read(src)?;
            let lower = html.to_ascii_lowercase();
            if let Some(head) = lower
                .find("<head")
                .and_then(|i| lower[i..].find('>').map(|j| i + j + 1))
            {
                html.insert_str(head, &meta(csp::RUN_HTML_CSP));
            }
            let lower = html.to_ascii_lowercase();
            if let Some(body) = lower
                .find("<body")
                .and_then(|i| lower[i..].find('>').map(|j| i + j + 1))
            {
                html.insert_str(body, &format!("\n<p role=\"note\">{warning}</p>\n"));
            }
            put(dst, &format!("<!-- {warning} -->\n{html}"))
        }
        _ => {
            if let Some(dir) = dst.parent() {
                io_at(dir, fs::create_dir_all(dir))?;
            }
            io_at(dst, fs::copy(src, dst)).map(|_| ())
        }
    }
}

fn fault(name: &str) -> bool {
    cfg!(feature = "fault-inject") && std::env::var("OG_LIVE_DEMO_FAULT").as_deref() == Ok(name)
}

fn frame_document(frames: Vec<Value>) -> Value {
    obj([
        ("interval_ms", Value::Int(INTERVAL_MS)),
        ("frames", Value::Array(frames)),
    ])
}

impl Demo<'_> {
    fn person(
        &self,
        out: &Path,
        run: &Run,
        frames: Vec<Value>,
        checked: Timestamp,
    ) -> Result<(), String> {
        let name = &run.name;
        let src = run.folder(&self.live.root);
        let dst = out.join("runs").join(name);
        let files = runs::run_files(name, Some(&src));
        let (base, events) = self.live.base_state(run);
        let state = self.live.base_json(&base, events.len());
        let mut results = plugins::results_body(&self.live, run, &results::FNS)?;
        let result = fields(&mut results)?;
        result.insert("generated_utc".into(), s(&checked.to_string()));
        if self.live.gated(name, &src, &state, None) {
            result.insert(
                "saved_attention".into(),
                self.live.attention(&src).0.unwrap_or(Value::Null),
            );
        }
        let html =
            pages::run_page(&self.live.cfg.assets, name, false).map_err(|e| e.to_string())?;
        put(&dst.join("index.html"), &page(html, "../../"))?;
        put_json(&dst.join("results.json"), &results)?;
        let mut extras = BTreeSet::new();
        for frame in &frames {
            if let Some(Value::Object(plugins)) =
                get(frame, "state").and_then(|st| get(st, "plugins"))
            {
                extras.extend(
                    plugins
                        .values()
                        .filter_map(|plugin| string(get(plugin, "file"))),
                );
            }
        }
        for rel in &files {
            let path = src.join(rel);
            if path.is_file() && inside(&path, &src) {
                copy_run_file(&path, &dst.join(rel), name)?;
            } else if !fault("omit-folder")
                && let Some(html) = runs::sole_html(run, &self.live.root, rel, &files)
            {
                copy_run_file(&path.join(html), &dst.join(rel).join("index.html"), name)?;
            }
        }
        for rel in extras {
            let path = src.join(rel);
            if !files.contains(rel) && path.is_file() && inside(&path, &src) {
                copy_run_file(&path, &dst.join(rel), name)?;
            }
        }
        put_json(&dst.join("frames.json"), &frame_document(frames))?;
        copy_run_file(
            &src.join("run_manifest.tsv"),
            &dst.join("run_manifest.tsv"),
            name,
        )?;
        copy_run_file(
            &src.join(format!("{name}_report.html")),
            &dst.join("report.html"),
            name,
        )?;
        for finding in items(get(&results, "findings")) {
            let id = string(get(finding, "id"))
                .filter(|id| plugins::valid_finding_id(id))
                .ok_or("a demo finding has an invalid id")?;
            let html =
                pages::finding_page(&self.live.cfg.assets, name, id).map_err(|e| e.to_string())?;
            put(
                &dst.join("finding").join(id).join("index.html"),
                &page(html, "../../../../"),
            )?;
        }
        put(&dst.join("viewer/index.html"), &page(
            "<!doctype html>\n<html lang=\"en\" data-theme=\"dark\"><head><meta charset=\"utf-8\">\
             <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\
             <title>Genome viewer: not in the demo</title><link rel=\"stylesheet\" href=\"/assets/live.css\"></head>\
             <body><header><h1>Genome viewer: not in the demo</h1></header><main>\
             <p>On your own live site this link opens the genome viewer on your own reads, on your own computer. \
             The demo has no reads to show.</p><p><a href=\"../\">Back to the run</a></p></main></body></html>\n".into(),
            "../../../",
        ))?;
        let input = dumps_compact(&obj([("results", results)]));
        let cache = self.live.cfg.state_dir.join("map");
        let doc = og(
            self.args,
            &[
                "map",
                "build",
                "--root",
                path_text(&self.live.root)?,
                "--run",
                name,
                "--assets",
                path_text(&self.live.cfg.repo_assets)?,
                "--cache-dir",
                path_text(&cache)?,
            ],
            input.as_bytes(),
        )?;
        io_at(&dst.join("map.json"), fs::write(dst.join("map.json"), doc))
    }

    fn site(&self, out: &Path, frames: Vec<Value>) -> Result<(), String> {
        let time = string(get(frames.last().ok_or("the demo has no frames")?, "t"))
            .ok_or("the last frame has no time")?;
        let checked: Timestamp = time.parse().map_err(|e: jiff::Error| e.to_string())?;
        let assets = &self.live.cfg.assets;
        let dst = out.join("assets");
        copy_tree(assets, &dst)?;
        for template in pages::PAGES {
            let path = dst.join(template);
            io_at(&path, fs::remove_file(&path))?;
        }
        let runs = self.live.runs();
        let listing = runs
            .values()
            .map(|r| (r.name.clone(), r.kind.clone()))
            .collect::<Vec<_>>();
        let home = pages::home_page(assets, &listing).map_err(|e| e.to_string())?
            .replace("=\"/runs/", "=\"runs/")
            .replace("=\"/runs.json\"", "=\"runs.json\"")
            .replace(
                "<form method=\"post\" action=\"/logout\" class=\"logout\"><button type=\"submit\">Log out</button></form>",
                "<p class=\"demo-note\">A demo: the people, the run and the results are invented.</p>",
            );
        put(&out.join("index.html"), &page(home, ""))?;
        let mut index = self.live.runs_index();
        if let Some(Value::Array(people)) = fields(&mut index)?.get_mut("runs") {
            for person in people {
                let person = fields(person)?;
                if let Some(end) = person.get("finished_utc").cloned() {
                    person.insert("last_heard_utc".into(), end);
                }
            }
        }
        relink(&mut index, "");
        put_json(&out.join("runs.json"), &index)?;
        self.person(out, &runs[SAMPLE], frames, checked)?;
        for (name, _, _) in OTHERS {
            let run = &runs[name];
            let one = frame(
                checked,
                self.live
                    .replay_state(run, checked.as_second() as f64, &DemoHost),
            );
            self.person(out, run, vec![one], checked)?;
        }

        let family = &runs[FAMILY];
        let family_src = family.folder(&self.live.root);
        let family_dst = out.join("runs").join(FAMILY);
        let mut family_json = plugins::family_page(&self.live, FAMILY, &family_src)?;
        relink(&mut family_json, "../../");
        put_json(&family_dst.join("family.json"), &family_json)?;
        let state = self
            .live
            .replay_state(family, checked.as_second() as f64, &DemoHost);
        put_json(
            &family_dst.join("frames.json"),
            &frame_document(vec![frame(checked, state)]),
        )?;
        let html = pages::run_page(assets, FAMILY, true).map_err(|e| e.to_string())?;
        put(&family_dst.join("index.html"), &page(html, "../../"))?;
        let people = family::members_of(&family_src)?;
        let files = out.join("family").join(FAMILY);
        put(
            &files.join("pedigree.svg"),
            &og_pedigree::layout::draw(&people, FAMILY)?,
        )?;
        for file in family::served_files(&family_src, &people) {
            copy_run_file(&family_src.join(&file), &files.join(&file), FAMILY)?;
        }

        let genes = out.join("map/genes");
        io_at(&genes, fs::create_dir_all(&genes))?;
        let cache = self.live.cfg.state_dir.join("map");
        let base = [
            "map",
            "genes",
            "--root",
            path_text(&self.live.root)?,
            "--cache-dir",
            path_text(&cache)?,
        ];
        for chrom in (1..=22)
            .map(|n| format!("chr{n}"))
            .chain(["chrX".into(), "chrY".into()])
        {
            let data = og(
                self.args,
                &[base.as_slice(), &["--chrom", &chrom]].concat(),
                b"",
            )?;
            let file = genes.join(format!("{chrom}.json"));
            io_at(&file, fs::write(&file, data))?;
        }
        let names = genes.join("names.json");
        let data = og(self.args, &[base.as_slice(), &["--names"]].concat(), b"")?;
        io_at(&names, fs::write(&names, data))
    }
}

fn build(args: &Args) -> Result<usize, String> {
    let mut events = Vec::new();
    let fixture = args.repo.join(FIXTURES).join("demo-001.weblog.jsonl");
    for line in read(&fixture)?.lines().filter(|l| !l.trim().is_empty()) {
        let event = ingest::parse(line.as_bytes())?;
        let time: Timestamp = string(get(&event, "utcTime"))
            .ok_or("weblog event has no utcTime")?
            .parse()
            .map_err(|e: jiff::Error| e.to_string())?;
        if events.last().is_some_and(|(previous, _)| previous > &time) {
            return Err("weblog events are out of order".into());
        }
        events.push((time, event));
    }
    let start = events.first().ok_or("the demo weblog is empty")?.0;
    let end = events.last().expect("nonempty weblog").0;
    if start.as_second() != RUN_STARTED {
        return Err(
            "the demo weblog start does not match the sample; run og dev demo-events".into(),
        );
    }

    let scratch;
    let work = match &args.work {
        Some(path) => {
            io_at(path, fs::create_dir_all(path))?;
            if io_at(path, fs::read_dir(path))?.next().is_some() {
                return Err(format!("--work {} is not empty", path.display()));
            }
            path
        }
        None => {
            scratch = Scratch::new("live-demo")?;
            &scratch.0
        }
    };
    let mut cfg = Config::new(work);
    cfg.repo_assets = args.repo.join("assets");
    cfg.steps = cfg.repo_assets.join("steps.tsv");
    cfg.assets = cfg.repo_assets.join("live");
    cfg.results = Some(results::FNS);
    cfg.ops.eta_cores = Some(16);
    cfg.ops.eta_box = Some("demo".into());
    cfg.og = args.og.clone();
    let steps = read_steps(&cfg.steps).map_err(|e| e.to_string())?;
    let log = io_at(&cfg.access_log, AccessLog::new(&cfg.access_log, 0.0))?;
    let mut live = Live::new(cfg, steps, log);
    live.tz = jiff::tz::TimeZone::UTC;
    let demo = Demo {
        args,
        live,
        start,
        end,
    };
    demo.sample(SAMPLE, "male", "1", args.records)?;
    for (name, sex, seed) in OTHERS {
        demo.sample(name, sex, seed, 2000.max(args.records / 4))?;
    }
    let family = demo.live.root.join(FAMILY);
    copy_tree(&demo.path("demo/DEMO-FAMILY"), &family)?;
    let couple = demo
        .live
        .root
        .join(SAMPLE)
        .join(format!("couple/{SAMPLE}_couple_carrier.tsv"));
    if couple.is_file() {
        io_at(
            &couple,
            fs::copy(&couple, family.join("couple_DEMO-001_DEMO-002.tsv")),
        )?;
    }
    demo.genes()?;
    let frames = demo.replay(&events)?;
    let count = frames.len();

    // Publish the site only after the full build succeeds. A sibling staging
    // directory keeps the rename on the output filesystem.
    let output = &args.output;
    let parent = output.parent().ok_or("--out has no parent folder")?;
    io_at(parent, fs::create_dir_all(parent))?;
    let staged = parent.join(format!(".og-live-demo-{}", std::process::id()));
    io_at(&staged, fs::create_dir(&staged))?;
    let staged = Scratch(staged);
    let site = staged.join("site");
    demo.site(&site, frames)?;
    if output.exists() {
        let empty = io_at(output, fs::read_dir(output))?.next().is_none();
        if !empty && !output.join("runs/DEMO-001/frames.json").is_file() {
            return Err(format!(
                "{} is not an empty folder or a built demo",
                output.display()
            ));
        }
        let previous = staged.join("previous");
        io_at(output, fs::rename(output, &previous))?;
        if let Err(error) = fs::rename(&site, output) {
            if let Err(restore) = fs::rename(&previous, output) {
                std::mem::forget(staged);
                return Err(format!(
                    "{}: {error}; could not restore {}: {restore}",
                    output.display(),
                    previous.display()
                ));
            }
            return Err(format!("{}: {error}", output.display()));
        }
    } else {
        io_at(output, fs::rename(&site, output))?;
    }
    Ok(count)
}

/// Write the complete static site, including the replay and the genome map.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if matches!(args.first().map(String::as_str), Some("-h" | "--help")) {
        let _ = writeln!(out, "{USAGE}");
        return 0;
    }
    let args = match Args::parse(args) {
        Ok(args) => args,
        Err(e) => {
            let _ = writeln!(err, "{USAGE}\nog live demo: {e}");
            return 2;
        }
    };
    match build(&args) {
        Ok(count) => {
            let _ = writeln!(
                out,
                "Demo written: {} ({count} frames)",
                args.output.display()
            );
            0
        }
        Err(e) => {
            let _ = writeln!(err, "og live demo: {e}");
            1
        }
    }
}
