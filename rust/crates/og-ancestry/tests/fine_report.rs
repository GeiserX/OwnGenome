//! `og ancestry-fine`, step 26b's commands, on made-up inputs, run as the
//! pipeline runs them (through the same entry point `og` calls):
//!
//! - `panel-map`: a group of fewer than 10 unrelated people (one Oceanian,
//!   three San) is no reference and no painting reference, and the install
//!   log names it as left out;
//! - `check-map`: panel-map's table passes; a short row, a flag other than
//!   0 or 1 and a table with no reference person each stop it;
//! - `chroms`: a chromosome is painted only with 100 SNVs spanning 5 cM;
//! - `report`: fractions with the bootstrap range, merged painting rows,
//!   the painted shares, a PNG of the right size, the left-out group.
//!
//! The unit tests of this crate and tests/parity check each command in
//! detail. Both pipeline reports of the table are checked by
//! `rust/crates/og-render/tests/ancestry_fine.rs`. These are the checks of
//! the retired `tests/test_ancestry_fine_report.py`, one for one.

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use flate2::read::ZlibDecoder;

/// A fresh folder for one test.
fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "og-ancestry-fine-report-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn p(path: &Path) -> String {
    path.display().to_string()
}

/// `og ancestry-fine ARGS`: (exit code, stdout, stderr).
fn og(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_ancestry::cli(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

fn rows(path: &Path) -> Vec<Vec<String>> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

const MAP_HEADER: &str = "sample\tpopulation\tcontinent\tregion\treference\tpaint\n";

#[test]
fn panel_map_and_check_map() {
    let d = scratch("panel-map");
    // 10 unrelated Yoruba, 3 unrelated San, one unrelated Papuan and two
    // related ones: in gnomAD's metadata 29 of the 30 Oceanians are related.
    let mut small: Vec<(String, &str, &str)> = Vec::new();
    small.extend((0..10).map(|i| (format!("y{i}"), "YRI", "false")));
    small.extend((0..3).map(|i| (format!("s{i}"), "San", "false")));
    small.push(("p0".into(), "PapuanHighlands", "false"));
    small.push(("p1".into(), "Bougainville", "true"));
    small.push(("p2".into(), "PapuanSepik", "true"));
    let mut meta = String::from("s\tpopulation\tsample_filters.all_samples_related\n");
    let mut bcf = String::new();
    for (s, pop, rel) in &small {
        meta += &format!("{s}\t{pop}\t{rel}\n");
        bcf += &format!("{s}\n");
    }
    fs::write(d.join("meta.tsv"), meta).unwrap();
    fs::write(d.join("bcf.txt"), bcf).unwrap();
    let samples = d.join("samples.tsv");
    let (code, out, err) = og(&[
        "panel-map",
        "--meta",
        &p(&d.join("meta.tsv")),
        "--bcf-samples",
        &p(&d.join("bcf.txt")),
        "--paint-per-group",
        "50",
        "--out",
        &p(&samples),
    ]);
    assert_eq!((code, err.as_str()), (0, ""), "panel-map: exits 0");
    assert!(
        out.contains(
            "Left out of the model, fewer than 10 unrelated reference people: \
             Oceania (continental and regional group), Southern Africa (regional group)"
        ),
        "panel-map: the install log names the groups left out: {out}"
    );
    let table = rows(&samples);
    let by: BTreeMap<&str, &Vec<String>> = table[1..].iter().map(|r| (r[0].as_str(), r)).collect();
    assert_eq!(
        table[0],
        [
            "sample",
            "population",
            "continent",
            "region",
            "reference",
            "paint"
        ],
        "panel-map: the header"
    );
    assert_eq!(
        by["p0"][2..],
        ["OCE", "OCE", "0", "0"],
        "panel-map: the one unrelated Oceanian keeps its groups but is no reference"
    );
    let san: std::collections::BTreeSet<&str> = (0..3)
        .map(|i| by[format!("s{i}").as_str()][4].as_str())
        .collect();
    assert_eq!(
        san.into_iter().collect::<Vec<_>>(),
        ["0"],
        "panel-map: a regional group of 3 (San) is no reference"
    );
    let west: std::collections::BTreeSet<(&str, &str)> = (0..10)
        .map(|i| {
            let r = by[format!("y{i}").as_str()];
            (r[4].as_str(), r[5].as_str())
        })
        .collect();
    assert_eq!(
        west.into_iter().collect::<Vec<_>>(),
        [("1", "1")],
        "panel-map: the 10 West Africans are references and painting references"
    );
    assert_eq!(
        og(&["check-map", "--map", &p(&samples)]).0,
        0,
        "check-map: panel-map's table passes"
    );
    for (name, text) in [
        ("short", format!("{MAP_HEADER}y0\tYRI\tAFR\n")),
        ("flag", format!("{MAP_HEADER}y10\tYRI\tAFR\tWAF\ttrue\t0\n")),
        ("empty", MAP_HEADER.to_owned()),
    ] {
        let f = d.join(format!("{name}.tsv"));
        fs::write(&f, text).unwrap();
        assert_eq!(
            og(&["check-map", "--map", &p(&f)]).0,
            1,
            "check-map: stops on the {name} table"
        );
    }
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn chroms_paints_only_with_100_snvs_over_5_cm() {
    let d = scratch("chroms");
    let maps = d.join("maps");
    fs::create_dir_all(&maps).unwrap();
    for c in ["chr20", "chr22", "chr19"] {
        // 1 cM per Mb
        fs::write(
            maps.join(format!("{c}.map")),
            format!("{c} . 0.0 1000000\n{c} . 10.0 11000000\n"),
        )
        .unwrap();
    }
    fs::write(
        d.join("counts.tsv"),
        "chr20\t500\t2000000\t8000000\nchr22\t500\t2000000\t4000000\n\
         chr21\t999\t1\t40000000\nchr19\t50\t1000000\t11000000\nhomref\t7\n",
    )
    .unwrap();
    let (_, out, _) = og(&[
        "chroms",
        "--counts",
        &p(&d.join("counts.tsv")),
        "--maps",
        &p(&maps),
    ]);
    assert_eq!(
        out, "chr20\n",
        "chroms: 6 cM painted; 2 cM, no map and 50 SNVs not"
    );
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn report_fractions_painting_png_and_text() {
    let d = scratch("report");
    let model = d.join("model");
    fs::create_dir_all(&model).unwrap();
    fs::write(
        model.join("levels.tsv"),
        "continent\t3\tEAS,WEA,AFR\nregion\t2\tIBE,ITA\n",
    )
    .unwrap();
    // The panel's samples.tsv, beside model/: an Oceanian the model left out.
    fs::write(
        d.join("samples.tsv"),
        format!(
            "{MAP_HEADER}i1\tIBS\tWEA\tIBE\t1\t1\nt1\tTSI\tWEA\tITA\t1\t1\n\
             p0\tPapuanHighlands\tOCE\tOCE\t0\t0\n"
        ),
    )
    .unwrap();
    let radm = d.join("admixture");
    fs::create_dir_all(&radm).unwrap();
    fs::write(radm.join("S_continent.3.Q"), "0.05 0.9 0.05\n").unwrap();
    let boot: String = (0..11)
        .map(|i| format!("0.05 0.{:02} 0.{:02}\n", 80 + i, 15 - i))
        .collect();
    fs::write(radm.join("S_continent.boot.Q"), boot).unwrap();
    fs::write(radm.join("S_region.2.Q"), "0.6 0.4\n").unwrap();
    let paint = d.join("paint");
    fs::create_dir_all(&paint).unwrap();
    fs::write(
        paint.join("chr20.msp.tsv"),
        "#Subpopulation order/codes: AFR=0\tEAS=1\tWEA=2\n\
         #chm\tspos\tepos\tsgpos\tegpos\tn snps\tS.0\tS.1\n\
         chr20\t1000\t5000\t1\t2\t10\t2\t2\nchr20\t5000\t9000\t2\t3\t10\t2\t0\n\
         chr20\t9000\t11000\t3\t4\t10\t2\t0\n",
    )
    .unwrap();
    // RFMix failed on chr19; chr20 was still painted.
    fs::write(paint.join("failed.txt"), "chr19\n").unwrap();
    let rep = d.join("out");
    let args = [
        "report".to_owned(),
        "--sample".into(),
        "S".into(),
        "--model".into(),
        p(&model),
        "--admixture".into(),
        p(&radm),
        "--paint".into(),
        p(&paint),
        "--panel".into(),
        "test".into(),
        "--out".into(),
        p(&rep),
    ];
    let argv: Vec<&str> = args.iter().map(String::as_str).collect();
    let (code, out, err) = og(&argv);
    assert_eq!((code, err.as_str()), (0, ""), "report: exits 0");
    let table = rows(&rep.join("S_ancestry_fine.tsv"));
    assert_eq!(
        table[0],
        ["level", "code", "group", "fraction", "se", "low", "high"],
        "report: the header"
    );
    let cont: Vec<&Vec<String>> = table.iter().filter(|r| r[0] == "continent").collect();
    assert_eq!(
        cont.iter().map(|r| r[1].as_str()).collect::<Vec<_>>(),
        ["WEA", "EAS", "AFR"],
        "report: continental groups, largest first"
    );
    assert_eq!(
        cont[0][3..],
        ["0.9000", "0.0332", "0.8025", "0.8975"],
        "report: West Eurasia's fraction and bootstrap range"
    );
    let region = table.iter().find(|r| r[0] == "region").unwrap();
    assert_eq!(
        region[4..],
        [".", ".", "."],
        "report: a level without bootstrap replicates has no range"
    );
    let bed = rows(&rep.join("S_painting.bed"));
    assert_eq!(
        bed[1..],
        [
            ["chr20", "999", "11000", "1", "WEA"],
            ["chr20", "999", "5000", "2", "WEA"],
            ["chr20", "4999", "11000", "2", "AFR"],
        ],
        "report: haplotype 1 one West Eurasian row, haplotype 2 split at 5000"
    );
    let shares: BTreeMap<&str, &str> = table
        .iter()
        .filter(|r| r[0] == "painting")
        .map(|r| (r[1].as_str(), r[3].as_str()))
        .collect();
    assert_eq!(
        shares,
        BTreeMap::from([("WEA", "0.7000"), ("AFR", "0.3000")]),
        "report: painted shares by length"
    );
    let png = fs::read(rep.join("S_painting.png")).unwrap();
    let w = u32::from_be_bytes(png[16..20].try_into().unwrap());
    let h = u32::from_be_bytes(png[20..24].try_into().unwrap());
    assert_eq!(
        (&png[..8], w, h),
        (&b"\x89PNG\r\n\x1a\n"[..], 1000, 6 + 22 * 19),
        "report: a PNG of 1000 pixels by the 22 chromosomes"
    );
    let find = |tag: &[u8]| png.windows(4).position(|x| x == tag).unwrap();
    let idat = &png[find(b"IDAT") + 4..find(b"IEND") - 8];
    let mut pixels = Vec::new();
    ZlibDecoder::new(idat).read_to_end(&mut pixels).unwrap();
    assert_eq!(
        pixels.len(),
        (h * (1 + 3 * w)) as usize,
        "report: the PNG's pixel data is whole"
    );
    let text = fs::read_to_string(rep.join("S_ancestry_fine.txt")).unwrap();
    assert_eq!(out, text, "report: stdout is the text");
    assert!(
        text.contains("West Eurasia (Europe, Middle East, North Africa)     90%  (80% to 90%)"),
        "report: the text names the group, the fraction and its range: {text}"
    );
    assert!(
        text.contains("not a nationality"),
        "report: the text says a name is ours, not a nationality: {text}"
    );
    assert!(
        text.contains(
            "Left out of the model: Oceania (continental and regional group). \
             The panel has fewer than 10 unrelated"
        ),
        "report: the text says which group the model left out and why: {text}"
    );

    fs::remove_file(paint.join("failed.txt")).unwrap();
    fs::remove_file(paint.join("chr20.msp.tsv")).unwrap();
    let (_, out, _) = og(&argv);
    assert!(
        out.contains("Chromosome painting: none."),
        "report: without a painting the text says none: {out}"
    );
    let _ = fs::remove_dir_all(&d);
}
