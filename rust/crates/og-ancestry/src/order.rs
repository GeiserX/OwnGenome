//! `order`: the group codes in the column order of the `.Q` (and `.P`)
//! file `admixture --supervised` wrote. Each column is the group whose
//! reference individuals have their largest mean share in it.

use crate::{Args, Fail, Result, Spec};
use std::io::Write;
use std::path::Path;

pub(crate) const SPEC: Spec = Spec {
    verb: "order",
    usage: "usage: og ancestry-fine order --pop FILE --q FILE\n",
    options: &["--pop", "--q"],
    required: &["--pop", "--q"],
};

/// The rows of numbers of a whitespace-separated matrix (`.Q`, `.P`),
/// blank lines skipped.
pub fn read_matrix(path: &Path) -> Result<Vec<Vec<f64>>> {
    let mut out = Vec::new();
    let mut n = 0u64;
    crate::for_each_line(path, |line| {
        n += 1;
        let mut row = Vec::new();
        for w in crate::split_ws(line) {
            let x = pycompat::num::parse_float(w).map_err(|_| {
                Fail(format!(
                    "{}: line {n}: not a number: {}",
                    path.display(),
                    crate::repr(w)
                ))
            })?;
            row.push(x);
        }
        if !row.is_empty() {
            out.push(row);
        }
        Ok(())
    })?;
    Ok(out)
}

/// The group of each column of a supervised ADMIXTURE `.Q`, given the
/// group of each of its rows (`-` or empty for a person in none).
pub fn column_order(pops: &[String], q: &[Vec<f64>]) -> std::result::Result<Vec<String>, String> {
    if pops.len() != q.len() {
        return Err(format!("{} .pop lines but {} .Q rows", pops.len(), q.len()));
    }
    let k = q.first().map_or(0, Vec::len);
    // A longer row's extra values are not read; a shorter one cannot be.
    if let Some(i) = q.iter().position(|r| r.len() < k) {
        return Err(format!(
            "row {} has {} values, the first {k}",
            i + 1,
            q[i].len()
        ));
    }
    let mut groups: Vec<&str> = pops
        .iter()
        .map(String::as_str)
        .filter(|p| *p != "-" && !p.is_empty())
        .collect();
    groups.sort_unstable();
    groups.dedup();
    if groups.len() != k {
        return Err(format!(
            "{} groups in the .pop file but K = {k}",
            groups.len()
        ));
    }
    let mut order: Vec<String> = Vec::with_capacity(k);
    for j in 0..k {
        // The first group (by name) with the largest mean wins a tie.
        let mut best: Option<(&str, f64)> = None;
        for &g in &groups {
            let vals = pops
                .iter()
                .zip(q)
                .filter(|(p, _)| *p == g)
                .map(|(_, r)| r[j]);
            let n = pops.iter().filter(|p| *p == g).count();
            let mean = pycompat::num::sum(vals) / n as f64;
            if best.is_none_or(|(_, m)| mean > m) {
                best = Some((g, mean));
            }
        }
        let (g, mean) = best.expect("k > 0 means at least one group");
        if mean < 0.5 {
            return Err(format!(
                "column {} has no group with a mean share above 0.5",
                j + 1
            ));
        }
        order.push(g.to_owned());
    }
    let mut seen = order.clone();
    seen.sort_unstable();
    seen.dedup();
    if seen.len() != k {
        return Err(format!(
            "two columns map to the same group: {}",
            crate::repr_list(&order)
        ));
    }
    Ok(order)
}

pub(crate) fn main(a: &Args, out: &mut dyn Write) -> Result<()> {
    let pops: Vec<String> = crate::read_lines(Path::new(a.req("--pop")))?
        .iter()
        .map(|l| pycompat::num::py_strip(l).to_owned())
        .collect();
    let q_path = a.req("--q");
    let q = read_matrix(Path::new(q_path))?;
    let order = column_order(&pops, &q).map_err(|e| Fail(format!("{q_path}: {e}")))?;
    writeln!(out, "{}", order.join(",")).map_err(|e| Fail(format!("cannot write stdout: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{Dir, run};

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| (*x).to_owned()).collect()
    }

    #[test]
    fn each_column_is_the_group_that_fills_it() {
        let pops = s(&["EUR", "EUR", "AFR", "EAS", "AFR"]);
        let q = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.99, 0.01, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.02, 0.98],
        ];
        assert_eq!(column_order(&pops, &q).unwrap(), ["EUR", "EAS", "AFR"]);
        let two: Vec<Vec<f64>> = q.iter().map(|r| r[..2].to_vec()).collect();
        assert_eq!(
            column_order(&pops, &two).unwrap_err(),
            "3 groups in the .pop file but K = 2"
        );
        let same = vec![vec![1.0, 0.0], vec![1.0, 0.0], vec![1.0, 0.0]];
        assert_eq!(
            column_order(&s(&["A", "A", "B"]), &same).unwrap_err(),
            "column 2 has no group with a mean share above 0.5"
        );
        let dup = vec![vec![0.9, 0.9], vec![0.1, 0.1]];
        assert_eq!(
            column_order(&s(&["A", "B"]), &dup).unwrap_err(),
            "two columns map to the same group: ['A', 'A']"
        );
    }

    #[test]
    fn the_command() {
        let d = Dir::new("order");
        let pop = d.put("ref.pop", "B\nA\n-\n");
        let q = d.put("ref.2.Q", "0.1 0.9\n\n0.8 0.2\n0.5 0.5\n");
        assert_eq!(
            run(&["order", "--pop", &pop, "--q", &q]),
            (0, "A,B\n".into(), String::new())
        );
        let bad = d.put("bad.Q", "0.1 x\n");
        let (rc, _, err) = run(&["order", "--pop", &pop, "--q", &bad]);
        assert_eq!(rc, 1);
        assert!(err.ends_with("bad.Q: line 1: not a number: 'x'\n"), "{err}");
    }
}
