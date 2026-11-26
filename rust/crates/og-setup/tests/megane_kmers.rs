//! `og setup megane-kmers` against MEGAnE's own `build_kmerset`.
//!
//! `tests/rust/vectors/megane-kmers/` holds small synthetic FASTAs, their
//! `.fai`, and the `.mk` and `.mi` MEGAnE's builder wrote for them in the
//! pinned `MEGANE_IMAGE` (its README has the recording steps). Each one is
//! built here under several thread counts and pass sizes and must come out
//! byte for byte the same. The FASTAs are what [`vector_fastas`] writes
//! (checked below), so they can be made again: `OG_WRITE_MEGANE_VECTORS=1
//! cargo test -p og-setup --test megane_kmers`.

use og_setup::megane_kmers::{self, Error, Genome, Plan};
use std::fs;
use std::path::{Path, PathBuf};

fn vectors() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../tests/rust/vectors/megane-kmers")
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-megane-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// xorshift64*: the same bases on every machine.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn seq(&mut self, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| b"ACGT"[(self.next() >> 62) as usize])
            .collect()
    }
}

fn rc(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'a' => b't',
            b't' => b'a',
            b'g' => b'c',
            b'c' => b'g',
            x => x,
        })
        .collect()
}

fn cat(parts: &[&[u8]]) -> Vec<u8> {
    parts.concat()
}

/// A FASTA of (name, sequence, line width); the file ends without a line
/// end when `last_newline` is false.
fn fasta(contigs: &[(&str, Vec<u8>, usize)], last_newline: bool) -> Vec<u8> {
    let mut out = Vec::new();
    for (name, seq, width) in contigs {
        out.extend_from_slice(format!(">{name} synthetic\n").as_bytes());
        for line in seq.chunks(*width) {
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    if !last_newline {
        out.pop();
    }
    out
}

/// One contig per input rule, each with the copies that make its k-mers
/// seen twice, so the rule shows in the output.
fn rules() -> Vec<u8> {
    let mut r = Rng(0x6d65_6761_6e65_0001);
    let mut c: Vec<(&str, Vec<u8>, usize)> = Vec::new();
    // Palindromes: P1 twice (one value per window, so kept), P2 once (not kept).
    let (h1, h2) = (r.seq(16), r.seq(16));
    let (p1, p2) = (cat(&[&h1, &rc(&h1)]), cat(&[&h2, &rc(&h2)]));
    let pal = cat(&[
        &r.seq(40),
        &p1,
        &r.seq(40),
        &p1,
        &r.seq(40),
        &p2,
        &r.seq(40),
    ]);
    c.push(("pal", pal, 60));
    // Both strands: S forward in one contig, its reverse complement in another.
    let s = r.seq(300);
    c.push(("fwd", s.clone(), 60));
    c.push(("rev", rc(&s), 50));
    // Lowercase is uppercase.
    let l = r.seq(120);
    c.push((
        "lower",
        cat(&[&r.seq(30), &l.to_ascii_lowercase(), &r.seq(30)]),
        60,
    ));
    c.push(("upper", cat(&[&r.seq(30), &l, &r.seq(30)]), 7));
    // N and n drop windows.
    let m = r.seq(150);
    let mut mn = m.clone();
    mn[70] = b'N';
    mn[100] = b'n';
    c.push(("nblock", mn, 60));
    c.push(("nfree", m, 60));
    // An N at base 0 alone keeps the first window, as A forward and as T
    // in the reverse complement; one at base 0 and base 5 does not.
    let q = r.seq(40);
    c.push(("n0", cat(&[b"N", &q]), 60));
    c.push(("a0", cat(&[&r.seq(20), b"A", &q, &r.seq(20)]), 60));
    c.push(("t0", cat(&[&r.seq(20), b"T", &q, &r.seq(20)]), 60));
    let q2 = r.seq(40);
    c.push(("n0lower", cat(&[b"n", &q2]), 60));
    c.push(("a0lower", cat(&[&r.seq(20), b"A", &q2, &r.seq(20)]), 60));
    let (q4, q3) = (r.seq(4), r.seq(40));
    c.push(("n0n5", cat(&[b"N", &q4, b"N", &q3]), 60));
    c.push((
        "a0a5",
        cat(&[&r.seq(20), b"A", &q4, b"A", &q3, &r.seq(20)]),
        60,
    ));
    // Other letters do not drop windows: forward they read as A, in the
    // reverse complement as T.
    let u = r.seq(260);
    let (mut ux, mut ua, mut ut) = (u.clone(), u.clone(), u);
    for (i, &x) in b"RYKMSWBDHVX-".iter().enumerate() {
        let p = 10 + i * 20;
        ux[p] = x;
        ua[p] = b'A';
        ut[p] = b'T';
    }
    c.push(("iupac", ux, 61));
    c.push(("iupac_a", ua, 60));
    c.push(("iupac_t", ut, 60));
    // No window spans two contigs: B's halves end one and start the next.
    let b = r.seq(32);
    c.push(("edge1", cat(&[&r.seq(40), &b[..16]]), 60));
    c.push(("edge2", cat(&[&b[16..], &r.seq(40)]), 60));
    c.push(("whole", cat(&[&r.seq(20), &b, &r.seq(20)]), 60));
    // Contigs of 31, 32 and 33 bases, each also inside a longer one.
    let (w1, w2, w3) = (r.seq(31), r.seq(32), r.seq(33));
    c.push(("s31", w1.clone(), 60));
    c.push(("s32", w2.clone(), 60));
    c.push(("s33", w3.clone(), 60));
    c.push((
        "shorts",
        cat(&[
            &r.seq(10),
            &w1,
            &r.seq(10),
            &w2,
            &r.seq(10),
            &w3,
            &r.seq(10),
        ]),
        60,
    ));
    // A lowercase n drops windows too: two copies of one contig with an n
    // at base 35, so a window over the n would be seen twice, and kept, if
    // the n did not drop it. Their own generator, so the contigs above and
    // below keep their bases.
    let mut nl = Rng(0x6d65_6761_6e65_0004).seq(70);
    nl[35] = b'n';
    c.push(("nlow1", nl.clone(), 60));
    c.push(("nlow2", nl, 60));
    // The last contig ends without a line end; its end is copied before it.
    let t = r.seq(137);
    c.push(("tailcopy", cat(&[&r.seq(10), &t[90..], &r.seq(10)]), 60));
    c.push(("tail", t, 60));
    fasta(&c, false)
}

/// A larger random genome with planted repeats: copies with a change every
/// 40 bases, reverse-complement copies, tandem repeats, poly-A, a
/// soft-masked stretch, N runs and other letters.
fn random() -> Vec<u8> {
    let mut r = Rng(0x6d65_6761_6e65_0002);
    let alu = r.seq(300);
    let mut chroms = Vec::new();
    for n in 0..4 {
        let mut s = Vec::new();
        for k in 0..12 {
            let gap = 400 + (r.next() % 400) as usize;
            s.extend(r.seq(gap));
            let mut copy = if (n + k) % 3 == 0 {
                rc(&alu)
            } else {
                alu.clone()
            };
            let start = (r.next() % 40) as usize;
            for i in (start..copy.len()).step_by(40) {
                copy[i] = b"ACGT"[(r.next() >> 62) as usize];
            }
            s.extend(copy);
            match k % 4 {
                0 => s.extend(b"CA".repeat(60)),
                1 => s.extend([b'A'; 80]),
                2 => s.extend([b'N'; 100]),
                _ => s.extend(r.seq(200).to_ascii_lowercase()),
            }
        }
        for i in (500..s.len()).step_by(997) {
            s[i] = b"RYKMN"[i % 5];
        }
        chroms.push(s);
    }
    let names = ["chrA", "chrB", "chrC", "chrD"];
    let c: Vec<(&str, Vec<u8>, usize)> = names
        .iter()
        .zip(chroms)
        .enumerate()
        .map(|(i, (n, s))| (*n, s, if i % 2 == 0 { 70 } else { 80 }))
        .collect();
    fasta(&c, true)
}

/// Two contigs; the hand-made `.fai` lists `y` twice.
fn twice() -> Vec<u8> {
    let mut r = Rng(0x6d65_6761_6e65_0003);
    fasta(&[("x", r.seq(100), 60), ("y", r.seq(80), 60)], true)
}

fn vector_fastas() -> Vec<(&'static str, Vec<u8>)> {
    vec![("rules", rules()), ("random", random()), ("twice", twice())]
}

#[test]
fn the_committed_fastas_are_what_the_generator_writes() {
    let dir = vectors();
    for (name, bytes) in vector_fastas() {
        let p = dir.join(format!("{name}.fa"));
        if std::env::var_os("OG_WRITE_MEGANE_VECTORS").is_some() {
            fs::write(&p, &bytes).unwrap();
        }
        assert!(
            fs::read(&p).unwrap() == bytes,
            "{} differs from the generator",
            p.display()
        );
    }
}

fn build(fa: &Path, out: &Path, plan: Plan) -> megane_kmers::Summary {
    let text = fs::read_to_string(format!("{}.fai", fa.display())).unwrap();
    let fai = megane_kmers::parse_fai(&text, "fai").unwrap();
    let g = Genome::load(fa, &fai).unwrap();
    megane_kmers::build_genome(&g, out, "og", plan).unwrap()
}

/// Every vector, under plans from one thread and one pass to many threads
/// and passes of a single value each.
#[test]
fn every_vector_is_byte_for_byte_megane_s_own_build() {
    let dir = vectors();
    for name in ["rules", "random", "twice"] {
        let mut plans = vec![(1, u64::MAX), (8, u64::MAX), (8, 64), (3, 500), (2, 4096)];
        let want_mk = fs::read(dir.join(format!("{name}.mk"))).unwrap();
        let want_mi = fs::read_to_string(dir.join(format!("{name}.mi"))).unwrap();
        assert!(
            !want_mk.is_empty(),
            "{name}.mk is empty: the vector shows nothing"
        );
        // One value a pass: every prefix gets a pass of its own, over the budget.
        if name != "random" {
            plans.push((1, 1));
        }
        for (threads, pass_values) in plans {
            let out = scratch(&format!("{name}-{threads}-{pass_values}"));
            let s = build(
                &dir.join(format!("{name}.fa")),
                &out,
                Plan {
                    threads,
                    pass_values,
                },
            );
            let got = fs::read(out.join("og.mk")).unwrap();
            assert!(
                got == want_mk,
                "{name} on {threads} threads, {pass_values} values a pass: .mk differs from MEGAnE's ({} values, MEGAnE {})",
                got.len() / 8,
                want_mk.len() / 8
            );
            assert_eq!(
                fs::read_to_string(out.join("og.mi")).unwrap(),
                want_mi,
                "{name} .mi"
            );
            assert_eq!(s.kept as usize, want_mk.len() / 8);
            fs::remove_dir_all(&out).unwrap();
        }
    }
}

/// The rules each vector must show, read off MEGAnE's own output, so a
/// vector that stopped exercising one fails here and not silently.
#[test]
fn the_vectors_exercise_every_rule() {
    let dir = vectors();
    let mk = fs::read(dir.join("rules.mk")).unwrap();
    let vals: Vec<u64> = mk
        .chunks_exact(8)
        .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
        .collect();
    assert!(vals.windows(2).all(|w| w[0] < w[1]), "strictly increasing");
    let code = |s: &[u8]| {
        s.iter().fold(0u64, |v, &b| {
            (v << 2)
                | match b {
                    b'T' => 1,
                    b'G' => 2,
                    b'C' => 3,
                    _ => 0,
                }
        })
    };
    let has = |s: &[u8]| vals.binary_search(&code(s)).is_ok();
    let fa = fs::read(dir.join("rules.fa")).unwrap();
    let seq = |name: &str| -> Vec<u8> {
        let head = format!(">{name} synthetic\n");
        let at = fa
            .windows(head.len())
            .position(|w| w == head.as_bytes())
            .unwrap()
            + head.len();
        let end = fa[at..]
            .iter()
            .position(|&b| b == b'>')
            .map_or(fa.len(), |e| at + e);
        fa[at..end]
            .iter()
            .copied()
            .filter(|&b| b != b'\n')
            .collect()
    };
    let pal = seq("pal");
    assert!(has(&pal[40..72]), "a palindrome seen twice is kept");
    assert!(!has(&pal[184..216]), "a palindrome seen once is not");
    let n0 = seq("n0");
    let a0: Vec<u8> = [b"A".as_slice(), &n0[1..32]].concat();
    assert!(has(&a0), "the N at base 0 reads as A forward");
    let t0: Vec<u8> = [b"T".as_slice(), &n0[1..32]].concat();
    assert!(has(&rc(&t0)), "and as T in the reverse complement");
    assert!(
        !has(&seq("a0a5")[20..52]),
        "an N at base 5 drops the first window"
    );
    let nl = seq("nlow1");
    assert_eq!(nl[35], b'n', "nlow1 has its n at base 35");
    assert!(has(&nl[0..32]), "a window before the lowercase n is kept");
    let over_n: Vec<u8> = [&nl[20..35], b"A".as_slice(), &nl[36..52]].concat();
    assert!(
        !has(&over_n) && !has(&rc(&[&nl[20..35], b"T".as_slice(), &nl[36..52]].concat())),
        "a lowercase n drops the windows that hold it"
    );
    let s = seq("fwd");
    let r = seq("rev");
    assert!(has(&s[0..32]) && has(&r[0..32]), "both strands");
    let u = seq("iupac_a");
    assert!(has(&u[0..32]), "a window over an IUPAC letter, read as A");
    assert!(
        has(&rc(&seq("iupac_t")[0..32])),
        "and as T in the reverse complement"
    );
    assert!(!has(&seq("whole")[20..52]), "no window spans two contigs");
    assert!(has(&seq("s32")), "a contig of 32 bases adds its one window");
    let t = seq("tail");
    assert!(
        has(&t[105..137]),
        "the last bases, after the last line without a line end"
    );
}

fn write(dir: &Path, name: &str, fa: &[u8], fai: &str) -> PathBuf {
    let p = dir.join(name);
    fs::write(&p, fa).unwrap();
    fs::write(format!("{}.fai", p.display()), fai).unwrap();
    p
}

fn refused(fa: &Path) -> String {
    let out = fa.parent().unwrap().join("out");
    match megane_kmers::build(fa, &out, "og", 1 << 30, 2) {
        Err(Error::Input(m)) => {
            assert!(!out.join("og.mk").exists() && !out.join("og.mi").exists());
            m
        }
        other => panic!("{} was not refused: {other:?}", fa.display()),
    }
}

#[test]
fn input_megane_reads_undefined_or_wrong_is_refused() {
    let d = scratch("refused");
    let seq = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT";
    let crlf = format!(">a\r\n{}\r\n{}\r\n", &seq[..20], &seq[20..]);
    let fa = write(&d, "crlf.fa", crlf.as_bytes(), "a\t36\t4\t20\t22\n");
    assert!(refused(&fa).contains("Windows (CRLF) line ends"));

    let mut high = format!(">a\n{seq}\n").into_bytes();
    high[10] = 0xC3;
    let fa = write(&d, "high.fa", &high, "a\t36\t3\t36\t37\n");
    assert!(refused(&fa).contains("byte 0xc3, outside ASCII"));

    let fa = write(
        &d,
        "short.fa",
        format!(">a\n{seq}\n").as_bytes(),
        "a\t60\t3\t36\t37\n",
    );
    assert!(refused(&fa).contains("does not describe"));

    let fa = write(
        &d,
        "stale.fa",
        format!(">a\n{}\n{}\n", &seq[..10], &seq[10..]).as_bytes(),
        "a\t36\t3\t20\t21\n",
    );
    assert!(refused(&fa).contains("is a line end"));

    let fa = write(&d, "gz.fa", b"\x1f\x8b\x08\x00rest", "a\t2\t3\t2\t3\n");
    assert!(refused(&fa).contains("gzip-compressed"));

    let fa = write(
        &d,
        "six.fa",
        format!(">a\n{seq}\n").as_bytes(),
        "a\t36\t3\t36\t37\t0\n",
    );
    assert!(refused(&fa).contains("a FASTA index line has 5"));

    let fa = write(
        &d,
        "tiny.fa",
        format!(">a\n{seq}\n").as_bytes(),
        "a\t36\t3\t36\t37\n",
    );
    let out = d.join("out");
    match megane_kmers::build(&fa, &out, "og", 1 << 20, 2) {
        Err(Error::Input(m)) => assert!(m.contains("--mem"), "{m}"),
        other => panic!("too little memory was not refused: {other:?}"),
    }
    // The message names the value as given, not rounded: --mem 0.25 on 8
    // threads needs 0.27 GB.
    match megane_kmers::build(&fa, &out, "og", 250_000_000, 8) {
        Err(Error::Input(m)) => assert!(
            m.starts_with("--mem 0.25 is too little") && m.ends_with("at least 0.3 GB"),
            "{m}"
        ),
        other => panic!("--mem 0.25 on 8 threads was not refused: {other:?}"),
    }
    fs::remove_dir_all(&d).unwrap();
}

/// A contig name that is not UTF-8: samtools indexes it and MEGAnE's
/// builder reads it, and og shows names only in its messages. `twice.fa`
/// with `x` renamed to the Latin-1 bytes `caf\xe9` must give `twice.mk`.
#[test]
fn a_contig_name_that_is_not_utf8_is_read() {
    let d = scratch("latin");
    let v = vectors();
    let fa = fs::read(v.join("twice.fa")).unwrap();
    let head = b">x synthetic\n";
    assert!(fa.starts_with(head), "twice.fa starts with contig x");
    let latin = [b">caf\xe9 synthetic\n".as_slice(), &fa[head.len()..]].concat();
    // The name is 3 bytes longer, so every offset moves by 3.
    let mut fai = Vec::new();
    for line in fs::read_to_string(v.join("twice.fa.fai")).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let name: &[u8] = if f[0] == "x" {
            b"caf\xe9"
        } else {
            f[0].as_bytes()
        };
        let offset: u64 = f[2].parse().unwrap();
        fai.extend_from_slice(name);
        fai.extend_from_slice(
            format!("\t{}\t{}\t{}\t{}\n", f[1], offset + 3, f[3], f[4]).as_bytes(),
        );
    }
    let p = d.join("latin.fa");
    fs::write(&p, &latin).unwrap();
    fs::write(d.join("latin.fa.fai"), &fai).unwrap();
    let out = d.join("out");
    let s = megane_kmers::build(&p, &out, "og", 1 << 30, 2).unwrap();
    assert!(
        fs::read(out.join("og.mk")).unwrap() == fs::read(v.join("twice.mk")).unwrap(),
        "the k-mers of twice.fa under another contig name"
    );
    assert_eq!(
        fs::read_to_string(out.join("og.mi")).unwrap(),
        fs::read_to_string(v.join("twice.mi")).unwrap()
    );
    assert_eq!(s.kept, 98);
    // A refusal still names the contig, with the byte shown as U+FFFD.
    let mut bad = latin.clone();
    bad[b">caf\xe9 synthetic\n".len() + 5] = 0xC3;
    fs::write(&p, &bad).unwrap();
    match megane_kmers::build(&p, &out, "og", 1 << 30, 2) {
        Err(Error::Input(m)) => assert!(
            m.contains("contig caf\u{FFFD}, base 6 is the byte 0xc3"),
            "{m}"
        ),
        other => panic!("a byte of 0x80 or more was not refused: {other:?}"),
    }
    fs::remove_dir_all(&d).unwrap();
}

#[test]
fn the_command_line_writes_both_files_and_keeps_old_ones_on_an_error() {
    let d = scratch("cli");
    let v = vectors();
    let fa = d.join("rules.fa");
    fs::copy(v.join("rules.fa"), &fa).unwrap();
    fs::copy(v.join("rules.fa.fai"), d.join("rules.fa.fai")).unwrap();
    let out = d.join("megane");
    let args = |extra: &[&str]| -> Vec<String> {
        let mut a: Vec<String> = ["--fa", fa.to_str().unwrap(), "--out", out.to_str().unwrap()]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        a.extend(extra.iter().map(|s| (*s).to_owned()));
        a
    };
    assert_eq!(
        megane_kmers::cli(&args(&["--threads", "2", "--mem", "1"])),
        0
    );
    // MEGAnE's default prefix is the FASTA's file name.
    assert_eq!(
        fs::read(out.join("rules.fa.mk")).unwrap(),
        fs::read(v.join("rules.mk")).unwrap()
    );
    assert_eq!(
        fs::read(out.join("rules.fa.mi")).unwrap(),
        fs::read(v.join("rules.mi")).unwrap()
    );
    assert_eq!(
        megane_kmers::cli(&args(&["--prefix", "p", "--threads", "0"])),
        2
    );
    assert_eq!(megane_kmers::cli(&args(&["--prefix", "a/b"])), 2);
    fs::write(&fa, b"\x1f\x8b").unwrap();
    assert_eq!(megane_kmers::cli(&args(&[])), 2);
    assert_eq!(
        fs::read(out.join("rules.fa.mk")).unwrap(),
        fs::read(v.join("rules.mk")).unwrap(),
        "a refused build leaves the old files"
    );
    assert_eq!(
        fs::read_dir(&out).unwrap().count(),
        2,
        "no temporary file left"
    );
    fs::remove_dir_all(&d).unwrap();
}
