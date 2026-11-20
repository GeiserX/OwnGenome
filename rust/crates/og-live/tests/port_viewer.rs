//! The local genome viewer (`og view serve`, `og_live::viewer`) on a
//! synthetic sample written here: a random 3 kb chr1 and 1 kb chr2
//! reference with its .fai, 30 reads on chr1 as a coordinate-sorted BAM with
//! its .bai, a bgzipped small-variant VCF with its .tbi and an SV VCF without
//! an index, each written byte by byte (no samtools or htslib).
//!
//! Moved here from tests/test_viewer.py, check for check: the config, the
//! locus check, the byte ranges, the refused paths and Host headers, the
//! loopback-only bind, and the page in a headless Chrome. New: a reversed
//! locus and a request a proxy forwarded are refused.
//!
//! The browser check needs Chrome and the igv.js release versions.env pins
//! (downloaded once into the test's scratch and checked by sha256). Without
//! Chrome it says SKIP, and fails in CI. Its red control serves an empty
//! igv.min.js: the page must then report an error, not "ready".

use flate2::Compression;
use flate2::Crc;
use flate2::write::DeflateEncoder;
use og_browser::{Browser, Scratch, Value, chrome_or_skip, repo_root};
use og_live::viewer::{self, Config, Logger, Server};
use pycompat::json;
use std::fs;
use std::io::{Read, Write};
use std::net::{TcpStream, UdpSocket};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex};
use std::time::Duration;

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut e = DeflateEncoder::new(Vec::new(), Compression::new(6));
    e.write_all(data).unwrap();
    let body = e.finish().unwrap();
    let mut crc = Crc::new();
    crc.update(data);
    let mut out = vec![0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff];
    out.extend(6u16.to_le_bytes());
    out.extend(b"BC");
    out.extend(2u16.to_le_bytes());
    out.extend(u16::try_from(body.len() + 25).unwrap().to_le_bytes());
    out.extend(&body);
    out.extend(crc.sum().to_le_bytes());
    out.extend(u32::try_from(data.len()).unwrap().to_le_bytes());
    out
}

fn reg2bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    for (shift, off) in [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> shift == end >> shift {
            return (off + (beg >> shift)) as u32;
        }
    }
    0
}

/// The per-reference part of a BAI or TBI: every record (ref, beg, end)
/// sits in one BGZF block that starts at `rec_voff`.
fn binning_index(
    n_ref: usize,
    recs: &[(usize, i64, i64)],
    rec_voff: i64,
    end_voff: i64,
) -> Vec<u8> {
    let mut out = Vec::new();
    for r in 0..n_ref {
        let mine: Vec<_> = recs.iter().filter(|x| x.0 == r).collect();
        if mine.is_empty() {
            out.extend(0i32.to_le_bytes());
            out.extend(0i32.to_le_bytes());
            continue;
        }
        let mut bins: Vec<u32> = mine.iter().map(|&&(_, b, e)| reg2bin(b, e)).collect();
        bins.sort_unstable();
        bins.dedup();
        out.extend((bins.len() as i32).to_le_bytes());
        for b in bins {
            out.extend(b.to_le_bytes());
            out.extend(1i32.to_le_bytes());
            out.extend(rec_voff.to_le_bytes());
            out.extend(end_voff.to_le_bytes());
        }
        let last = mine.iter().map(|x| x.2).max().unwrap();
        let n_intv = ((last - 1) >> 14) + 1;
        out.extend((n_intv as i32).to_le_bytes());
        for _ in 0..n_intv {
            out.extend(rec_voff.to_le_bytes());
        }
    }
    out
}

/// A BAM of `reads` (ref index, 0-based pos, sequence), sorted, 50M each,
/// with its .bai.
fn write_bam(path: &Path, contigs: &[(String, String)], reads: &[(usize, i64, String)]) {
    let mut text = "@HD\tVN:1.6\tSO:coordinate\n".to_owned();
    for (n, s) in contigs {
        text.push_str(&format!("@SQ\tSN:{n}\tLN:{}\n", s.len()));
    }
    text.push_str("@RG\tID:S\tSM:S\n");
    let mut head = b"BAM\x01".to_vec();
    head.extend((text.len() as i32).to_le_bytes());
    head.extend(text.as_bytes());
    head.extend((contigs.len() as i32).to_le_bytes());
    for (n, s) in contigs {
        head.extend((n.len() as i32 + 1).to_le_bytes());
        head.extend(n.as_bytes());
        head.push(0);
        head.extend((s.len() as i32).to_le_bytes());
    }
    let code = |c: u8| "=ACMGRSVTWYHKDBN".bytes().position(|x| x == c).unwrap() as u8;
    let (mut body, mut recs) = (Vec::new(), Vec::new());
    for (i, (r, pos, seq)) in reads.iter().enumerate() {
        let name = format!("r{i}\0");
        let sb = seq.as_bytes();
        let packed: Vec<u8> = sb
            .chunks(2)
            .map(|c| (code(c[0]) << 4) | c.get(1).map_or(0, |&x| code(x)))
            .collect();
        let len = sb.len() as i64;
        let mut rec = Vec::new();
        rec.extend((*r as i32).to_le_bytes());
        rec.extend((*pos as i32).to_le_bytes());
        rec.push(name.len() as u8);
        rec.push(60);
        rec.extend((reg2bin(*pos, pos + len) as u16).to_le_bytes());
        rec.extend(1u16.to_le_bytes());
        rec.extend(0u16.to_le_bytes());
        rec.extend((len as i32).to_le_bytes());
        for v in [-1i32, -1, 0] {
            rec.extend(v.to_le_bytes());
        }
        rec.extend(name.as_bytes());
        rec.extend(((len as u32) << 4).to_le_bytes());
        rec.extend(&packed);
        rec.extend(std::iter::repeat_n(30u8, sb.len()));
        rec.extend(b"RGZS\0");
        body.extend((rec.len() as i32).to_le_bytes());
        body.extend(rec);
        recs.push((*r, *pos, pos + len));
    }
    let (b0, b1) = (bgzf_block(&head), bgzf_block(&body));
    let mut all = b0.clone();
    all.extend(&b1);
    all.extend(BGZF_EOF);
    fs::write(path, all).unwrap();
    let mut bai = b"BAI\x01".to_vec();
    bai.extend((contigs.len() as i32).to_le_bytes());
    bai.extend(binning_index(
        contigs.len(),
        &recs,
        (b0.len() as i64) << 16,
        ((b0.len() + b1.len()) as i64) << 16,
    ));
    fs::write(format!("{}.bai", path.display()), bai).unwrap();
}

/// A bgzipped VCF of `records` (chrom, 1-based pos, ref, alt, info), with
/// a .tbi when `index`.
fn write_vcf(
    path: &Path,
    contigs: &[(String, String)],
    records: &[(&str, i64, String, String, &str)],
    index: bool,
) {
    let mut head = "##fileformat=VCFv4.2\n".to_owned();
    for (n, s) in contigs {
        head.push_str(&format!("##contig=<ID={n},length={}>\n", s.len()));
    }
    head.push_str("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n");
    head.push_str("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n");
    head.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    head.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n");
    let body: String = records
        .iter()
        .map(|(c, p, r, a, i)| format!("{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t{i}\tGT\t0/1\n"))
        .collect();
    let (b0, b1) = (bgzf_block(head.as_bytes()), bgzf_block(body.as_bytes()));
    let mut all = b0.clone();
    all.extend(&b1);
    all.extend(BGZF_EOF);
    fs::write(path, all).unwrap();
    if index {
        let order: Vec<&str> = contigs.iter().map(|(n, _)| n.as_str()).collect();
        let mut names: Vec<&str> = records.iter().map(|r| r.0).collect();
        names.sort_by_key(|n| order.iter().position(|o| o == n));
        names.dedup();
        let nm: Vec<u8> = names.iter().flat_map(|n| n.bytes().chain([0])).collect();
        let recs: Vec<_> = records
            .iter()
            .map(|(c, p, r, _, _)| {
                let i = names.iter().position(|n| n == c).unwrap();
                (i, p - 1, p - 1 + r.len() as i64)
            })
            .collect();
        let mut tbi = b"TBI\x01".to_vec();
        for v in [
            names.len() as i32,
            2,
            1,
            2,
            0,
            i32::from(b'#'),
            0,
            nm.len() as i32,
        ] {
            tbi.extend(v.to_le_bytes());
        }
        tbi.extend(&nm);
        tbi.extend(binning_index(
            names.len(),
            &recs,
            (b0.len() as i64) << 16,
            ((b0.len() + b1.len()) as i64) << 16,
        ));
        let mut out = bgzf_block(&tbi);
        out.extend(BGZF_EOF);
        fs::write(format!("{}.tbi", path.display()), out).unwrap();
    }
}

/// A random sequence, the same on every run.
fn sequence(seed: &mut u64, n: usize) -> String {
    (0..n)
        .map(|_| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(*seed >> 62) as usize] as char
        })
        .collect()
}

/// The genome folder `gd`: reference/synthetic.fasta and sample S. The
/// FASTA path and the contigs.
fn make_sample(gd: &Path) -> (PathBuf, Vec<(String, String)>) {
    let mut seed = 7;
    let contigs = vec![
        ("chr1".to_owned(), sequence(&mut seed, 3000)),
        ("chr2".to_owned(), sequence(&mut seed, 1000)),
    ];
    fs::create_dir_all(gd.join("reference")).unwrap();
    let fa = gd.join("reference/synthetic.fasta");
    let (mut fasta, mut fai, mut off) = (String::new(), String::new(), 0);
    for (n, s) in &contigs {
        off += n.len() + 2;
        fai.push_str(&format!("{n}\t{}\t{off}\t60\t61\n", s.len()));
        let lines: Vec<&str> = s
            .as_bytes()
            .chunks(60)
            .map(|c| std::str::from_utf8(c).unwrap())
            .collect();
        fasta.push_str(&format!(">{n}\n{}\n", lines.join("\n")));
        off += s.len() + lines.len();
    }
    fs::write(&fa, fasta).unwrap();
    fs::write(gd.join("reference/synthetic.fasta.fai"), fai).unwrap();
    let s = gd.join("S");
    for d in ["aligned", "vcf", "sv_merged"] {
        fs::create_dir_all(s.join(d)).unwrap();
    }
    let seq = &contigs[0].1;
    let reads: Vec<_> = (400..700)
        .step_by(10)
        .map(|p| (0, p as i64, seq[p..p + 50].to_owned()))
        .collect();
    write_bam(&s.join("aligned/S_sorted.bam"), &contigs, &reads);
    let r = &seq[450..451];
    let alt = "ACGT".replace(r, "")[..1].to_owned();
    write_vcf(
        &s.join("vcf/S.vcf.gz"),
        &contigs,
        &[("chr1", 451, r.to_owned(), alt, ".")],
        true,
    );
    write_vcf(
        &s.join("sv_merged/S_sv_consensus.vcf.gz"),
        &contigs,
        &[(
            "chr1",
            1000,
            seq[999..1000].to_owned(),
            "<DEL>".into(),
            "SVTYPE=DEL;END=1500",
        )],
        false,
    );
    (fa, contigs)
}

fn sample_config(gd: &Path, fa: &Path, locus: Option<&str>) -> Config {
    let s = gd.join("S");
    viewer::build_config(
        "S",
        fa,
        Some(&s.join("aligned/S_sorted.bam")),
        Some(&s.join("vcf/S.vcf.gz")),
        Some(&s.join("sv_merged/S_sv_consensus.vcf.gz")),
        locus,
    )
    .unwrap()
}

/// The answer to one request: status, headers (lower-case names) and body.
fn get(
    port: u16,
    method: &str,
    path: &str,
    headers: &[(&str, &str)],
) -> (u16, Vec<(String, String)>, Vec<u8>) {
    let mut c = TcpStream::connect(("127.0.0.1", port)).unwrap();
    c.set_read_timeout(Some(Duration::from_secs(20))).unwrap();
    let mut req = format!("{method} {path} HTTP/1.1\r\n");
    if !headers.iter().any(|(k, _)| k.eq_ignore_ascii_case("host")) {
        req.push_str(&format!("Host: 127.0.0.1:{port}\r\n"));
    }
    for (k, v) in headers {
        req.push_str(&format!("{k}: {v}\r\n"));
    }
    req.push_str("\r\n");
    c.write_all(req.as_bytes()).unwrap();
    let mut raw = Vec::new();
    c.read_to_end(&mut raw).unwrap();
    let end = raw
        .windows(4)
        .position(|w| w == b"\r\n\r\n")
        .expect("a response head");
    let head = String::from_utf8(raw[..end].to_vec()).unwrap();
    let mut lines = head.split("\r\n");
    let code = lines
        .next()
        .unwrap()
        .split(' ')
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    let hs = lines
        .filter_map(|l| l.split_once(':'))
        .map(|(k, v)| (k.trim().to_ascii_lowercase(), v.trim().to_owned()))
        .collect();
    (code, hs, raw[end + 4..].to_vec())
}

fn header<'a>(h: &'a [(String, String)], name: &str) -> Option<&'a str> {
    h.iter().find(|(k, _)| k == name).map(|(_, v)| v.as_str())
}

fn assets() -> PathBuf {
    repo_root().join("assets/viewer")
}

/// A stand-in igv.min.js for the HTTP checks, which never run it.
fn fake_igv(dir: &Path) -> PathBuf {
    let p = dir.join("igv.min.js");
    fs::write(&p, "// igv stand-in\n").unwrap();
    p
}

fn file_of(cfg: &Config, url: &str) -> PathBuf {
    cfg.files[url.split('/').nth(1).unwrap().parse::<usize>().unwrap()].clone()
}

fn str_at<'a>(v: &'a Value, path: &[&str]) -> Option<&'a str> {
    let mut v = v;
    for k in path {
        v = match v {
            Value::Object(m) => m.get(*k)?,
            Value::Array(a) => a.get(k.parse::<usize>().ok()?)?,
            _ => return None,
        };
    }
    match v {
        Value::Str(s) => Some(s),
        _ => None,
    }
}

fn tracks(v: &Value) -> Vec<Value> {
    match v {
        Value::Object(m) => match m.get("tracks") {
            Some(Value::Array(a)) => a.clone(),
            _ => Vec::new(),
        },
        _ => Vec::new(),
    }
}

#[test]
fn the_config_names_the_samples_files_and_checks_the_locus() {
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "viewer-config");
    let gd = scratch.path().join("genome");
    let (fa, _) = make_sample(&gd);
    let cfg = sample_config(&gd, &fa, Some("chr1:451"));
    assert_eq!(viewer::BIND, "127.0.0.1");
    assert_eq!(cfg.locus, "chr1:451");
    assert!(
        cfg.files[0].ends_with("reference/synthetic.fasta")
            && cfg.files[1].ends_with("synthetic.fasta.fai")
    );
    assert_eq!(
        str_at(&cfg.json, &["reference", "fastaURL"]),
        Some("data/0/synthetic.fasta")
    );
    let t = tracks(&cfg.json);
    let kinds: Vec<(String, String, String)> = t
        .iter()
        .map(|x| {
            (
                str_at(x, &["type"]).unwrap().to_owned(),
                str_at(x, &["format"]).unwrap().to_owned(),
                file_of(&cfg, str_at(x, &["url"]).unwrap())
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .into_owned(),
            )
        })
        .collect();
    assert_eq!(
        kinds,
        [
            ("alignment".into(), "bam".into(), "S_sorted.bam".into()),
            ("variant".into(), "vcf".into(), "S.vcf.gz".into()),
            (
                "variant".into(),
                "vcf".into(),
                "S_sv_consensus.vcf.gz".into()
            ),
        ]
    );
    assert!(file_of(&cfg, str_at(&t[0], &["indexURL"]).unwrap()).ends_with("S_sorted.bam.bai"));
    assert!(
        str_at(&t[1], &["indexURL"])
            .unwrap()
            .ends_with("S.vcf.gz.tbi")
    );
    let Value::Object(sv) = &t[2] else { panic!() };
    assert_eq!(
        sv.get("indexed"),
        Some(&Value::Bool(false)),
        "the SV VCF without an index loads whole"
    );
    // The same bytes Python's json.dumps(cfg, indent=1) wrote.
    assert!(cfg.to_json().starts_with("{\n \"sample\": \"S\",\n \"locus\": \"chr1:451\",\n \"reference\": {\n  \"id\": \"reference\""));
    assert_eq!(json::loads(&cfg.to_json()).unwrap(), cfg.json);

    assert_eq!(
        sample_config(&gd, &fa, None).locus,
        "chr1",
        "no locus opens the first contig"
    );
    for bad in [
        "chr9:100",
        "chr1:99999",
        "chr1;alert(1)",
        "BRCA1 x",
        "chr1:2000-1000",
    ] {
        assert!(
            viewer::build_config("S", &fa, None, None, None, Some(bad)).is_err(),
            "locus {bad} was accepted"
        );
    }
    let bam = gd.join("S/aligned/S_sorted.bam");
    fs::rename(
        gd.join("S/aligned/S_sorted.bam.bai"),
        scratch.path().join("bai"),
    )
    .unwrap();
    let e = viewer::build_config("S", &fa, Some(&bam), None, None, None).unwrap_err();
    assert!(e.contains("samtools index"), "a BAM without an index: {e}");
    let e = viewer::build_config(
        "S",
        &gd.join("reference/none.fasta"),
        None,
        None,
        None,
        None,
    )
    .unwrap_err();
    assert!(e.contains("needs its FASTA and .fai"), "{e}");
}

#[test]
fn manta_diploid_svs_come_before_the_consensus() {
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "viewer-files");
    let s = scratch.path().join("S");
    fs::create_dir_all(s.join("sv_merged")).unwrap();
    fs::write(s.join("sv_merged/S_sv_consensus.vcf.gz"), "").unwrap();
    let f = viewer::sample_files(&s, "S");
    assert_eq!((f.reads, f.vcf), (None, None));
    assert_eq!(
        f.sv_vcf,
        Some(s.join("sv_merged/S_sv_consensus.vcf.gz")),
        "the consensus when Manta's is missing"
    );
    fs::create_dir_all(s.join("manta/results/variants")).unwrap();
    fs::write(s.join("manta/results/variants/diploidSV.vcf.gz"), "").unwrap();
    fs::create_dir_all(s.join("aligned")).unwrap();
    fs::write(s.join("aligned/S_sorted.cram"), "").unwrap();
    let f = viewer::sample_files(&s, "S");
    assert_eq!(
        f.sv_vcf,
        Some(s.join("manta/results/variants/diploidSV.vcf.gz"))
    );
    assert_eq!(f.reads, Some(s.join("aligned/S_sorted.cram")));
    fs::write(s.join("aligned/S_sorted.bam"), "").unwrap();
    assert_eq!(
        viewer::sample_files(&s, "S").reads,
        Some(s.join("aligned/S_sorted.bam")),
        "the BAM before the CRAM"
    );
}

/// An address of this machine that is not loopback, if it has one.
fn other_address() -> Option<String> {
    let so = UdpSocket::bind("0.0.0.0:0").ok()?;
    so.connect("192.0.2.1:9").ok()?;
    let ip = so.local_addr().ok()?.ip().to_string();
    (!ip.starts_with("127.") && ip != "0.0.0.0").then_some(ip)
}

#[test]
fn the_server_answers_only_what_the_config_names() {
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "viewer-http");
    let gd = scratch.path().join("genome");
    let (fa, _) = make_sample(&gd);
    let cfg = sample_config(&gd, &fa, Some("chr1:451"));
    let srv = Server::start(&cfg, &assets(), &fake_igv(scratch.path()), 0, None).unwrap();
    let port = srv.port();

    let (code, h, body) = get(port, "GET", "/", &[]);
    let csp = header(&h, "content-security-policy").unwrap_or_default();
    assert!(
        code == 200 && String::from_utf8_lossy(&body).contains("igv.min.js"),
        "{code}"
    );
    assert!(
        csp.contains("connect-src 'self'")
            && csp.contains("default-src 'none'")
            && csp.contains("frame-ancestors 'none'"),
        "{csp}"
    );
    assert_eq!(csp, viewer::CSP);
    assert_eq!(header(&h, "x-content-type-options"), Some("nosniff"));
    assert_eq!(header(&h, "cache-control"), Some("no-store"));

    let (code, h, body) = get(port, "GET", "/config.json", &[]);
    assert_eq!(
        (code, header(&h, "content-type")),
        (200, Some("application/json"))
    );
    let sc = json::loads(std::str::from_utf8(&body).unwrap()).unwrap();
    assert_eq!(sc, cfg.json);
    let st = tracks(&sc);
    let names: Vec<&str> = st
        .iter()
        .map(|t| str_at(t, &["url"]).unwrap().rsplit('/').next().unwrap())
        .collect();
    assert_eq!(names, ["S_sorted.bam", "S.vcf.gz", "S_sv_consensus.vcf.gz"]);

    let raw = fs::read(gd.join("S/aligned/S_sorted.bam")).unwrap();
    let bam = format!("/{}", str_at(&tracks(&sc)[0], &["url"]).unwrap());
    let (code, h, body) = get(port, "GET", &bam, &[("Range", "bytes=10-109")]);
    assert_eq!((code, body.as_slice()), (206, &raw[10..110]));
    assert_eq!(
        header(&h, "content-range").map(str::to_owned),
        Some(format!("bytes 10-109/{}", raw.len()))
    );
    let (code, _, body) = get(port, "GET", &bam, &[("Range", "bytes=-28")]);
    assert_eq!(
        (code, body.as_slice()),
        (206, &BGZF_EOF[..]),
        "a suffix range (the BGZF EOF block)"
    );
    let (code, h, body) = get(port, "GET", &bam, &[]);
    assert_eq!(
        (code, body.len(), header(&h, "accept-ranges")),
        (200, raw.len(), Some("bytes"))
    );
    let past = format!("bytes={}-", raw.len() + 5);
    let (code, h, _) = get(port, "GET", &bam, &[("Range", &past)]);
    assert_eq!(code, 416);
    assert_eq!(
        header(&h, "content-range").map(str::to_owned),
        Some(format!("bytes */{}", raw.len()))
    );
    let (code, h, body) = get(port, "HEAD", &bam, &[("Range", "bytes=0-9")]);
    assert_eq!(
        (code, header(&h, "content-length"), body.len()),
        (206, Some("10"), 0),
        "HEAD sends no body"
    );

    for p in [
        "/data/0/../../reference/synthetic.fasta.fai",
        "/data/2/..%2F..%2Fvcf%2FS.vcf.gz",
        "/data/99/S_sorted.bam",
        "/serve.py",
        "/data/2/S.vcf.gz",
    ] {
        assert_eq!(get(port, "GET", p, &[]).0, 404, "{p} must be 404");
    }
    assert_eq!(
        get(
            port,
            "GET",
            "/config.json",
            &[("Host", "attacker.example:8765")]
        )
        .0,
        403,
        "a foreign Host"
    );
    assert_eq!(
        get(port, "GET", "/config.json", &[("Host", "localhost:8765")]).0,
        200
    );
    assert_eq!(get(port, "POST", "/config.json", &[]).0, 405);
    // D12: the viewer has no password, so a request a proxy forwarded is refused.
    for (k, v) in [
        ("X-Forwarded-For", "203.0.113.9"),
        ("Forwarded", "for=203.0.113.9"),
        ("X-Forwarded-Host", "genome.example"),
        ("X-Forwarded-Proto", "https"),
        ("X-Real-IP", "203.0.113.9"),
    ] {
        let (code, _, body) = get(port, "GET", "/config.json", &[(k, v)]);
        assert_eq!(code, 403, "a request with {k} must be refused");
        assert!(
            String::from_utf8_lossy(&body).contains("ssh"),
            "the refusal says how to reach it"
        );
    }

    if let Some(ip) = other_address() {
        assert!(
            TcpStream::connect_timeout(
                &format!("{ip}:{port}").parse().unwrap(),
                Duration::from_secs(2)
            )
            .is_err(),
            "the viewer must not be reachable on {ip}"
        );
    } else {
        println!("[SKIP] no non-loopback address to test the bind against");
    }
    drop(srv);

    // A file name with a space, '#' or '?' (a REF_FASTA the user named so).
    let odd = scratch.path().join("odd");
    fs::create_dir_all(&odd).unwrap();
    fs::copy(&fa, odd.join("my ref#1?.fasta")).unwrap();
    fs::copy(
        gd.join("reference/synthetic.fasta.fai"),
        odd.join("my ref#1?.fasta.fai"),
    )
    .unwrap();
    let ocfg =
        viewer::build_config("S", &odd.join("my ref#1?.fasta"), None, None, None, None).unwrap();
    let url = str_at(&ocfg.json, &["reference", "fastaURL"])
        .unwrap()
        .to_owned();
    assert_eq!(url, "data/0/my%20ref%231%3F.fasta");
    let osrv = Server::start(&ocfg, &assets(), &fake_igv(scratch.path()), 0, None).unwrap();
    let (code, _, body) = get(
        osrv.port(),
        "GET",
        &format!("/{url}"),
        &[("Range", "bytes=0-4")],
    );
    assert_eq!((code, body.as_slice()), (206, &b">chr1"[..]));
}

#[test]
fn og_view_serve_refuses_a_bad_start() {
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "viewer-cli");
    let gd = scratch.path().join("genome");
    let (fa, _) = make_sample(&gd);
    let igv = fake_igv(scratch.path());
    // Every run names a port that is taken, so a check that stops failing
    // ends with exit 1 instead of a server that never returns.
    let busy = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
    let busy_port = busy.local_addr().unwrap().port().to_string();
    let run = |extra: &[&str]| {
        let mut args: Vec<String> = ["--sample", "S", "--fasta"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        args.push(fa.display().to_string());
        args.extend([
            "--igv-js".into(),
            igv.display().to_string(),
            "--assets".into(),
            assets().display().to_string(),
        ]);
        args.extend(["--port".into(), busy_port.clone()]);
        args.extend(extra.iter().map(|s| (*s).to_owned()));
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = viewer::cli(&args, &mut out, &mut err);
        (code, String::from_utf8(err).unwrap())
    };
    let (code, err) = run(&["--locus", "chr1:2000-1000"]);
    assert_eq!(code, 2);
    assert!(
        err.starts_with("ERROR: 'chr1:2000-1000' ends before it starts"),
        "{err}"
    );
    let (code, err) = run(&["--port", "http"]);
    assert!(code == 2 && err.contains("is not a port number"), "{err}");
    let (code, err) = run(&["--og-scaffold-probe"]);
    assert!(code == 2 && err.contains("unrecognized arguments"), "{err}");
    let (code, err) = run(&["--assets", "/nonexistent"]);
    assert!(
        code == 2 && err.contains("the viewer page not found"),
        "{err}"
    );
    // A port that is taken: exit 1 with the variable to set.
    let (code, err) = run(&[]);
    assert!(
        code == 1 && err.contains("Set OWNGENOME_VIEW_PORT"),
        "{err}"
    );
}

/// The igv.js release versions.env pins, downloaded once into `dir` and
/// checked by sha256.
fn pinned_igv(dir: &Path) -> PathBuf {
    let env = fs::read_to_string(repo_root().join("versions.env")).unwrap();
    let pin = |k: &str| {
        env.lines()
            .find_map(|l| l.strip_prefix(&format!("{k}=\"")))
            .and_then(|v| v.strip_suffix('"'))
            .unwrap_or_else(|| panic!("versions.env sets no {k}"))
            .to_owned()
    };
    let (version, sha) = (pin("IGV_JS_VERSION"), pin("IGV_JS_SHA256"));
    let path = dir.join(format!("igv-{version}.min.js"));
    if !path.is_file() {
        let url = format!("https://cdn.jsdelivr.net/npm/igv@{version}/dist/igv.min.js");
        let ok = Command::new("curl")
            .args(["-fsSL", "--retry", "3", "-o"])
            .arg(&path)
            .arg(&url)
            .status()
            .is_ok_and(|s| s.success());
        assert!(ok, "cannot download {url}");
    }
    let sum = Command::new("shasum")
        .args(["-a", "256"])
        .arg(&path)
        .output()
        .or_else(|_| Command::new("sha256sum").arg(&path).output())
        .expect("shasum or sha256sum");
    let got = String::from_utf8_lossy(&sum.stdout)
        .split_whitespace()
        .next()
        .unwrap_or_default()
        .to_owned();
    if got != sha {
        let _ = fs::remove_file(&path);
        panic!("igv.min.js has sha256 {got}, versions.env pins {sha}");
    }
    path
}

/// Opens the viewer of `cfg` in Chrome with `igv` as igv.min.js: the
/// page's status, what the policy blocked, and the requests it made.
fn load_page(
    chrome: &Path,
    scratch: &Path,
    name: &str,
    cfg: &Config,
    igv: &Path,
) -> (String, String, Vec<String>) {
    let lines = Arc::new(Mutex::new(Vec::new()));
    let sink = Arc::clone(&lines);
    let log: Logger = Box::new(move |l: &str| sink.lock().unwrap().push(l.to_owned()));
    let srv = Server::start(cfg, &assets(), igv, 0, Some(log)).unwrap();
    let mut b = Browser::launch(chrome, &scratch.join(format!("chrome-{name}"))).unwrap();
    b.go(&format!("http://127.0.0.1:{}/#chr1:451", srv.port()))
        .unwrap();
    let done = "(function(){var s=document.documentElement.dataset.pgpStatus||'';return s==='ready'||s.indexOf('error')===0;})()";
    b.wait(done, Duration::from_secs(60));
    // igv.js reads the reads after it reports ready; give it the time.
    b.wait(
        "performance.getEntriesByType('resource').some(function(e){return /S_sorted\\.bam$/.test(e.name);})",
        Duration::from_secs(20),
    );
    let text = |v: Result<Value, String>| match v {
        Ok(Value::Str(s)) => s,
        _ => String::new(),
    };
    let status = text(b.eval("document.documentElement.dataset.pgpStatus || ''"));
    let blocked = text(b.eval("document.documentElement.dataset.pgpBlocked || ''"));
    drop(b);
    drop(srv);
    let log = lines.lock().unwrap().clone();
    (status, blocked, log)
}

#[test]
fn the_page_starts_igv_on_the_config() {
    let Some(chrome) = chrome_or_skip("load the viewer page") else {
        return;
    };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "viewer-browser");
    let gd = scratch.path().join("genome");
    let (fa, _) = make_sample(&gd);
    let cfg = sample_config(&gd, &fa, Some("chr1:451"));
    let igv = pinned_igv(Path::new(env!("CARGO_TARGET_TMPDIR")));

    let (status, blocked, log) = load_page(&chrome, scratch.path(), "real", &cfg, &igv);
    let all = log.join("\n");
    assert_eq!(
        status, "ready",
        "igv.js did not start on the config; requests:\n{all}"
    );
    assert_eq!(blocked, "", "the policy blocked what the page needed");
    assert!(
        all.contains("S_sorted.bam.bai")
            && log
                .iter()
                .any(|l| l.starts_with("\"GET /data/") && l.contains("/S_sorted.bam HTTP")),
        "igv.js read the BAM's index and a range of the BAM:\n{all}"
    );
    assert!(
        all.contains("synthetic.fasta") && all.contains("S.vcf.gz"),
        "igv.js read the reference and the VCF:\n{all}"
    );
    assert!(
        !log.iter().any(|l| l.contains("\" 4")),
        "every request of the page was served:\n{all}"
    );

    // Red control: an empty igv.min.js. The same check must see an error.
    let empty = scratch.path().join("empty-igv.min.js");
    fs::write(&empty, "").unwrap();
    let (status, _, _) = load_page(&chrome, scratch.path(), "control", &cfg, &empty);
    assert!(
        status.starts_with("error: "),
        "control: the page reported {status:?} without igv.js"
    );
}
