use super::*;

fn lines(text: &str) -> Vec<String> {
    text.split_inclusive('\n').map(str::to_owned).collect()
}

fn members(text: &str) -> (Members, Vec<String>) {
    let (rows, errors) = parse_pedigree("p.tsv", &lines(text));
    assert!(errors.is_empty(), "{errors:?}");
    build_members(&rows, &HashMap::new(), "p.tsv", "the samplesheet")
}

/// The twin of `test_kinship` in tests/test_family_report.py.
#[test]
fn expected_relatedness() {
    let (m, errors) = members(
        "sample\tfather\tmother\tpartner\n\
         gp\t.\t.\t.\ngm\t.\t.\t.\n\
         f\tgp\tgm\t.\nu\tgp\tgm\t.\nm\t.\t.\tf\nw\t.\t.\t.\nw2\t.\t.\t.\n\
         a\tf\tm\t.\nb\tf\tm\t.\nh\tf\tw\t.\nc\tu\tw2\t.\n",
    );
    assert!(errors.is_empty(), "{errors:?}");
    let mut kin = Kinship::new(&m);
    for (x, y, want, label) in [
        ("a", "f", 0.5, "parent-child"),
        ("a", "b", 0.5, "full siblings"),
        ("a", "h", 0.25, "half siblings"),
        ("a", "gp", 0.25, "second degree"),
        ("a", "c", 0.125, "third degree"),
        ("f", "m", 0.0, "partners"),
        ("m", "w", 0.0, "unrelated"),
    ] {
        let got = 2.0 * kin.phi(x, y);
        assert!((got - want).abs() < 1e-9, "{x} {y}: {got}");
        assert_eq!(declared_relationship(x, y, &m, got), label, "{x} {y}");
    }
}

#[test]
fn pedigree_problems_are_named() {
    let (_, e) = parse_pedigree("p.tsv", &lines("kid\tdad\n"));
    assert!(
        e[0].contains("must be a header with a 'sample' column"),
        "{e:?}"
    );
    let (_, e) = parse_pedigree("p.tsv", &lines("sample\tcolour\n"));
    assert_eq!(
        e[0],
        "p.tsv: unknown column colour; the columns are family, sample, sex, father, mother, partner"
    );
    let (_, e) = parse_pedigree("p.tsv", &[]);
    assert_eq!(e, ["p.tsv is empty"]);

    let all = |text: &str| {
        let (rows, _) = parse_pedigree("p.tsv", &lines(text));
        build_members(&rows, &HashMap::new(), "p.tsv", "the samplesheet").1
    };
    let e = all("sample\tsex\tfather\nkid\tmale\tmum\nmum\tfemale\t.\n");
    assert_eq!(e, ["'mum' is the father of 'kid' but its sex is female"]);
    let e = all("sample\tfather\nkid\tkid\n");
    assert_eq!(e, ["p.tsv line 2: 'kid' is named as its own father"]);
    let e = all("sample\tfather\na\tb\nb\ta\n");
    assert!(
        e.iter().all(|l| l.contains("among its own ancestors")) && e.len() == 2,
        "{e:?}"
    );
    let e = all("sample\nkid\nkid\n");
    assert_eq!(
        e,
        ["p.tsv line 3: sample 'kid' has a second row; each person has one"]
    );

    let (rows, _) = parse_pedigree("p.tsv", &lines("sample\tsex\nkid\tmale\n"));
    let declared: HashMap<String, String> = [("kid".to_owned(), "female".to_owned())].into();
    let (_, e) = build_members(&rows, &declared, "p.tsv", "the samplesheet");
    assert_eq!(
        e,
        ["p.tsv line 2: the pedigree says 'kid' is male, the samplesheet says female"]
    );
}

#[test]
fn unnamed_relatives_get_a_member() {
    let (m, errors) = members("sample\tsex\tfather\tmother\nkid\tmale\tdad\tmum\n");
    assert!(errors.is_empty());
    assert_eq!(m.keys().collect::<Vec<_>>(), ["kid", "dad", "mum"]);
    assert_eq!(m["dad"].sex.as_deref(), Some("male"));
    assert!(!m["dad"].row && m["kid"].row);
}

#[test]
fn families_split_by_column() {
    let (rows, _) = parse_pedigree(
        "p.tsv",
        &lines("family\tsample\tfather\nfa\tkid\tdad\nfa\tdad\t.\nfb\tx\t.\n"),
    );
    let (f, e) = split_families(&rows, "p.tsv");
    assert!(e.is_empty(), "{e:?}");
    assert_eq!(
        f.keys().cloned().collect::<Vec<_>>(),
        [Some("fa".to_owned()), Some("fb".to_owned())]
    );
    assert_eq!(f[&Some("fa".to_owned())].len(), 2);
    let (rows, _) = parse_pedigree("p.tsv", &lines("family\tsample\nfa\tkid\n\tdad\n"));
    let (_, e) = split_families(&rows, "p.tsv");
    assert!(!e.is_empty());
}

#[test]
fn par_is_bed() {
    let dir = std::env::temp_dir().join(format!("og-pedigree-par-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("par.bed");
    std::fs::write(&p, "chrX\t10\t20\nshort\n").unwrap();
    let par = Par::read(p.to_str()).unwrap();
    assert!(!par.contains("chrX", 10) && par.contains("chrX", 11) && par.contains("chrX", 20));
    assert!(!par.contains("chrX", 21) && !par.contains("chrY", 15));
    assert!(Par::read(None).unwrap().is_empty());
    std::fs::write(&p, "chrX\tten\t20\n").unwrap();
    assert!(Par::read(p.to_str()).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}
