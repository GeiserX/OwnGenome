//! Throughput of og-io's readers and writer and of og-vcf's record view on a
//! VCF, and a byte-for-byte check of everything they read and write.
//!
//!     cargo run --release -p og-vcf --example throughput -- FILE.vcf.gz [WORKERS...]
//!
//! Loads the file once as the reference, then streams it again through each
//! reader and compares every byte with the reference as it goes, so the
//! timings hold no allocation. Prints one Markdown table row per measure (MB/s
//! of uncompressed VCF) and exits 1 when any reader, the round trip or the
//! writer gives other bytes.

use std::io::{self, BufRead, Read, Write};
use std::num::NonZero;
use std::process::{Command, ExitCode, Stdio};
use std::time::Instant;

use og_io::bgzf;
use og_vcf::{Header, Lines, Record, write_info};

fn mbps(bytes: usize, secs: f64) -> String {
    format!("{:.0}", bytes as f64 / secs / 1e6)
}

fn verdict(same: bool) -> &'static str {
    if same { "same bytes" } else { "DIFFERENT" }
}

fn timed<T>(f: impl FnOnce() -> io::Result<T>) -> io::Result<(T, f64)> {
    let t = Instant::now();
    let out = f()?;
    Ok((out, t.elapsed().as_secs_f64()))
}

/// Compares what is written to it with `want`, without keeping it.
struct Check<'a> {
    want: &'a [u8],
    pos: usize,
    same: bool,
}

impl<'a> Check<'a> {
    fn new(want: &'a [u8]) -> Self {
        Check {
            want,
            pos: 0,
            same: true,
        }
    }

    fn take(&mut self, got: &[u8]) {
        let end = self.pos + got.len();
        self.same &= self.want.get(self.pos..end) == Some(got);
        self.pos = end;
    }

    fn done(&self) -> bool {
        self.same && self.pos == self.want.len()
    }
}

/// Streams a reader into a [`Check`].
fn stream(mut r: impl BufRead, want: &[u8]) -> io::Result<bool> {
    let mut check = Check::new(want);
    loop {
        let buf = r.fill_buf()?;
        if buf.is_empty() {
            return Ok(check.done());
        }
        let n = buf.len();
        check.take(buf);
        r.consume(n);
    }
}

/// Parses every record, rebuilds INFO from its items and writes the line
/// back into a reused buffer that is checked against `want`.
fn round_trip(input: impl BufRead, want: &[u8]) -> io::Result<(usize, bool)> {
    let mut check = Check::new(want);
    let mut lines = Lines::new(input);
    let mut out = Vec::with_capacity(1 << 16);
    Header::read(&mut lines)?.write(&mut out, &[])?;
    check.take(&out);
    let mut info = Vec::new();
    let mut records = 0;
    while let Some(line) = lines.next_line()? {
        let rec = Record::parse(line).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        info.clear();
        write_info(rec.info_items(), &mut info);
        out.clear();
        rec.write_with_info(&info, &mut out);
        out.push(b'\n');
        check.take(&out);
        records += 1;
    }
    Ok((records, check.done()))
}

fn run() -> io::Result<bool> {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: throughput FILE.vcf.gz [WORKERS...]");
    let mut counts: Vec<usize> = args
        .map(|a| a.parse().expect("WORKERS is a number"))
        .collect();
    if counts.is_empty() {
        counts = vec![1, 2, 4, bgzf::workers().get()];
        counts.dedup();
    }
    let compressed = std::fs::metadata(&path)?.len();
    let mut want = Vec::new();
    og_io::open(&path)?.read_to_end(&mut want)?;
    let size = want.len();
    let mut ok = true;

    println!("| measure | MB/s | check |");
    println!("|---|---|---|");
    println!(
        "| input: {} ({:.0} MB gzip, {:.0} MB VCF), {} cores | | |",
        path.rsplit('/').next().unwrap_or(&path),
        compressed as f64 / 1e6,
        size as f64 / 1e6,
        bgzf::workers()
    );

    let (same, secs) = timed(|| stream(og_io::open(&path)?, &want))?;
    ok &= same;
    println!(
        "| og_io::open, one thread (flate2) | {} | {} |",
        mbps(size, secs),
        verdict(same)
    );

    if let Ok((same, secs)) = timed(|| {
        let mut child = Command::new("gzip")
            .arg("-dc")
            .arg(&path)
            .stdout(Stdio::piped())
            .spawn()?;
        let same = stream(
            io::BufReader::with_capacity(1 << 16, child.stdout.take().expect("piped")),
            &want,
        )?;
        Ok(child.wait()?.success() && same)
    }) {
        ok &= same;
        println!("| gzip -dc | {} | {} |", mbps(size, secs), verdict(same));
    }

    for &n in &counts {
        let w = NonZero::new(n).expect("WORKERS above 0");
        let (same, secs) = timed(|| stream(bgzf::open(&path, w)?, &want))?;
        ok &= same;
        println!(
            "| og_io::bgzf::open, {n} workers | {} | {} |",
            mbps(size, secs),
            verdict(same)
        );
    }

    let ((records, same), secs) = timed(|| round_trip(io::Cursor::new(&want[..]), &want))?;
    ok &= same;
    println!(
        "| og_vcf parse + INFO rebuild + write, from memory, one thread ({records} records) | {} | {} |",
        mbps(size, secs),
        verdict(same)
    );

    let n = *counts.last().expect("a worker count");
    let w = NonZero::new(n).expect("WORKERS above 0");
    let ((records, same), secs) = timed(|| round_trip(bgzf::open(&path, w)?, &want))?;
    ok &= same;
    println!(
        "| bgzf::open ({n} workers) + og_vcf round trip ({records} records) | {} | {} |",
        mbps(size, secs),
        verdict(same)
    );

    for &n in &[1, n] {
        let w = NonZero::new(n).expect("WORKERS above 0");
        let (file, secs) = timed(|| {
            let mut wr = bgzf::Writer::new(Vec::with_capacity(size / 8), w)?;
            wr.write_all(&want)?;
            wr.finish()
        })?;
        let out_mb = file.len() as f64 / 1e6;
        let same = stream(og_io::reader(io::Cursor::new(file))?, &want)?;
        ok &= same;
        println!(
            "| bgzf::Writer, {n} workers ({out_mb:.0} MB out) | {} | {} |",
            mbps(size, secs),
            if same {
                "reads back the same"
            } else {
                "DIFFERENT"
            }
        );
    }
    Ok(ok)
}

fn main() -> ExitCode {
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => {
            eprintln!("ERROR: some output differs from the input");
            ExitCode::FAILURE
        }
        Err(e) => {
            eprintln!("ERROR: {e}");
            ExitCode::FAILURE
        }
    }
}
