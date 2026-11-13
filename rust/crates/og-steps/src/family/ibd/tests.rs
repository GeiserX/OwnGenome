use super::*;

fn map(rows: &[(i64, f64)]) -> GeneticMap {
    GeneticMap::from_rows("t.map", rows.to_vec()).unwrap()
}

#[test]
fn dosage_reads_two_biallelic_alleles_only() {
    assert_eq!(dosage("0|0"), 0);
    assert_eq!(dosage("0/1:35"), 1);
    assert_eq!(dosage("1|0"), 1);
    assert_eq!(dosage("1|1"), 2);
    for bad in ["./.", "0|2", "1", "0||1", "0/1/1", "", ".|1", "01"] {
        assert_eq!(dosage(bad), -1, "{bad}");
    }
}

#[test]
fn bisect_is_pythons_on_sorted_and_unsorted_lists() {
    let a = [1, 3, 3, 5];
    assert_eq!((bisect_left(&a, 3), bisect_right(&a, 3)), (1, 3));
    assert_eq!((bisect_left(&a, 0), bisect_right(&a, 9)), (0, 4));
    // Python's own answers on an unsorted list.
    let u = [5, 1, 4, 2, 3];
    assert_eq!((bisect_left(&u, 3), bisect_right(&u, 3)), (2, 2));
}

#[test]
fn map_interpolates_and_holds_flat_outside() {
    let m = map(&[(300, 3.0), (100, 1.0), (200, 1.0)]);
    assert_eq!(m.at(50), 1.0);
    assert_eq!(m.at(150), 1.0);
    assert_eq!(m.at(250), 2.0);
    assert_eq!(m.at(400), 3.0);
    assert!(GeneticMap::from_rows("t.map", vec![(1, 0.0)]).is_err());
}

#[test]
fn union_and_intersect() {
    assert_eq!(
        union(vec![(5, 9), (1, 3), (3, 4), (8, 12)]),
        [(1, 4), (5, 12)]
    );
    assert_eq!(
        intersect(&[(1, 10), (20, 30)], &[(5, 25)]),
        [(5, 10), (20, 25)]
    );
    assert!(intersect(&[(1, 5)], &[(5, 9)]).is_empty());
}

#[test]
fn ibd2_needs_both_haplotypes_of_a_shared() {
    // a1-b1 over 0-100 and a2-b2 over 50-150: both copies over 50-100.
    assert_eq!(ibd2(&[(1, 1, 0, 100), (2, 2, 50, 150)]), [(50, 100)]);
    // a1-b1 and a1-b2: one haplotype of a only, never both copies.
    assert!(ibd2(&[(1, 1, 0, 100), (1, 2, 50, 150)]).is_empty());
    assert_eq!(ibd2(&[(1, 2, 0, 100), (2, 1, 60, 90)]), [(60, 90)]);
}

/// Ten markers at 10, 20, ... 100; `opp` lists the ones where the two are
/// opposite homozygotes.
fn pair(opp: &[i64]) -> (Vec<i64>, Vec<Dosage>, Vec<Dosage>) {
    let pos: Vec<i64> = (1..=10).map(|i| i * 10).collect();
    let dx = vec![0; 10];
    let dy = pos
        .iter()
        .map(|p| if opp.contains(p) { 2 } else { 0 })
        .collect();
    (pos, dx, dy)
}

#[test]
fn stitch_joins_across_one_opposite_homozygote_not_two() {
    let (pos, dx, dy) = pair(&[50]);
    assert_eq!(stitch(&[(20, 40), (60, 90)], &pos, &dx, &dy), [(10, 100)]);
    let (pos2, dx2, dy2) = {
        // Two opposite homozygotes between the stretches.
        let pos: Vec<i64> = vec![10, 20, 30, 40, 45, 50, 60, 70, 80, 90, 100];
        let dy: Vec<Dosage> = pos
            .iter()
            .map(|p| if [45, 50].contains(p) { 2 } else { 0 })
            .collect();
        (pos, vec![0; 11], dy)
    };
    assert_eq!(
        stitch(&[(20, 40), (60, 90)], &pos2, &dx2, &dy2),
        [(10, 40), (60, 100)],
        "two opposite homozygotes keep the stretches apart"
    );
    // Reaching the ends: one opposite homozygote before, none after.
    let (pos3, dx3, dy3) = pair(&[10]);
    assert_eq!(stitch(&[(30, 80)], &pos3, &dx3, &dy3), [(10, 100)]);
}

#[test]
fn stitch_does_not_reach_an_end_across_two() {
    let (pos, dx, dy) = pair(&[10, 20]);
    assert_eq!(stitch(&[(30, 80)], &pos, &dx, &dy), [(30, 100)]);
    assert!(stitch(&[], &pos, &dx, &dy).is_empty());
}

#[test]
fn degree_tells_siblings_from_parent_and_child() {
    assert_eq!(degree(0.5, 0.0, 1), "first degree (parent-child)");
    assert_eq!(degree(0.5, 0.09, 1), "first degree (parent-child)");
    assert_eq!(degree(0.625, 0.25, 3), "first degree (full siblings)");
    assert_eq!(degree(0.5, 0.10, 3), "first degree (full siblings)");
    assert_eq!(degree(0.8, 0.9, 1), "same person or identical twins");
    assert_eq!(degree(0.05, 0.0, 1), "fourth degree");
    assert_eq!(
        degree(0.01, 0.0, 1),
        "distant or unrelated (short segments only)"
    );
    assert_eq!(degree(0.0, 0.0, 0), "no shared segment");
}

#[test]
fn join_is_os_path_join_for_a_file_name() {
    assert_eq!(join("maps", "chr1.map"), "maps/chr1.map");
    assert_eq!(join("maps/", "chr1.map"), "maps/chr1.map");
    assert_eq!(join("", "chr1.map"), "chr1.map");
}

/// A scratch folder under the target dir's temp, removed on drop.
struct Scratch(std::path::PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!("og-ibd-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }
    fn path(&self, p: &str) -> String {
        self.0.join(p).display().to_string()
    }
    fn write(&self, p: &str, text: &str) {
        let f = self.0.join(p);
        std::fs::create_dir_all(f.parent().unwrap()).unwrap();
        std::fs::write(f, text).unwrap();
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn args(d: &Scratch) -> Args {
    let v: Vec<String> = [
        "--pedigree",
        &d.path("pedigree.tsv"),
        "--chroms",
        &d.path("chroms.txt"),
        "--maps",
        &d.path("maps"),
        "--phased",
        &d.path("phased"),
        "--ibd",
        &d.path("ibd"),
        "--out-dir",
        &d.path("out"),
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    super::super::args::parse(&SPEC, &v).ok().unwrap()
}

/// A trio-sized family on one chromosome: A and B share it all.
fn family(d: &Scratch) {
    d.write(
        "pedigree.tsv",
        "sample\tsex\tfather\tmother\tpartner\nA\tmale\tB\t.\t.\n",
    );
    d.write("chroms.txt", "chr1\n");
    d.write(
        "maps/chr1.map",
        "chr1\t.\t0.0\t1\nchr1\t.\t10.0\t10000000\n",
    );
    let mut vcf = String::from(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
    );
    for i in 1..10 {
        vcf.push_str(&format!(
            "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n",
            i * 1_000_000
        ));
    }
    d.write("phased/chr1.vcf.gz", &vcf);
    d.write(
        "ibd/chr1.ibd.gz",
        "B\t1\tA\t2\tchr1\t2000000\t8000000\t6.0\n",
    );
}

fn stop(r: Result<()>) -> Die {
    match r {
        Err(Stop::Die(d)) => d,
        other => panic!("expected a stop, got {other:?}"),
    }
}

use super::super::Stop;

#[test]
fn a_parent_and_child_share_the_whole_chromosome() {
    let d = Scratch::new("pc");
    family(&d);
    run(&args(&d)).unwrap();
    let pairs = std::fs::read_to_string(d.path("out/ibd_pairs.tsv")).unwrap();
    assert_eq!(
        pairs.lines().nth(1).unwrap(),
        "A\tB\tparent-child\t0.500\t1\t8.00\t0.00\t8.00\t8.00\t1.0000\t0.500\tfirst degree (parent-child)"
    );
    let segs = std::fs::read_to_string(d.path("out/ibd_segments.tsv")).unwrap();
    assert_eq!(
        segs.lines().nth(1).unwrap(),
        "A\tB\tchr1\t1000000\t9000000\t8.00\t0.00"
    );
}

#[test]
fn a_missing_input_is_an_error_line_and_exit_2() {
    for gone in [
        "maps/chr1.map",
        "phased/chr1.vcf.gz",
        "ibd/chr1.ibd.gz",
        "chroms.txt",
        "pedigree.tsv",
    ] {
        let d = Scratch::new("gone");
        family(&d);
        std::fs::remove_file(d.path(gone)).unwrap();
        let die = stop(run(&args(&d)));
        assert_eq!(die.code, 2, "{gone}");
        assert!(
            die.lines[0].starts_with(&format!("cannot read {}", d.path(gone))),
            "{gone}: {:?}",
            die.lines
        );
        assert!(
            !Path::new(&d.path("out/ibd_pairs.tsv")).exists(),
            "{gone}: no output"
        );
    }
}

#[test]
fn a_short_record_or_a_listed_twice_chromosome_stops() {
    let d = Scratch::new("short");
    family(&d);
    d.write(
        "phased/chr1.vcf.gz",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\nchr1\t5\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n",
    );
    let die = stop(run(&args(&d)));
    assert!(
        die.lines[0].ends_with("line 2: 1 sample columns, the #CHROM line names 2"),
        "{:?}",
        die.lines
    );
    family(&d);
    d.write("chroms.txt", "chr1\n chr1 \n");
    let die = stop(run(&args(&d)));
    assert!(
        die.lines[0].ends_with("line 2: chr1 is listed twice"),
        "{:?}",
        die.lines
    );
}

#[test]
fn a_phased_record_of_another_chromosome_stops() {
    let d = Scratch::new("otherchrom");
    family(&d);
    d.write(
        "phased/chr1.vcf.gz",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\nchr1\t1000000\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\nchr2\t5000000\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n",
    );
    let die = stop(run(&args(&d)));
    assert!(
        die.lines[0].ends_with("phased/chr1.vcf.gz line 3: chromosome chr2 in the file of chr1"),
        "{:?}",
        die.lines
    );
}
