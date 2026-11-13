#!/usr/bin/env bash
# summary-collect-full.sh: writes the sample folders of the whole-file
# summary-collect parity cases (og summary collect without --only) in the
# current directory, the same every time. It runs the generators of the
# three groups (summary-collect.sh, summary-clinical.sh and
# summary-collect-c.sh) and joins their folders, so one folder holds every
# section and the post pass reads across the groups: the repeat readings
# with the CPIC state, Locityper's primary calls with pypgx and Cyrius,
# CNVpytor's chrY filter with step 33's declared sex, the cyp2d6 block,
# not_run and not_assessed. On top: a run with failed steps, an empty
# folder, previous runs kept in history/ (one half-copied, one whose
# summary is no JSON) with a database re-check, and a telomere table whose
# infinite read length stopped the Python collector (a golden: og keeps
# the summary and marks that one section unreadable).
# py/ holds the modules the oracle imports at the oracle's commit, with the
# real run_delta.py for recheck and new_since. The case stages
# cb_vep.vcf.gz and its index beside the folders.
# scripts/ci/oracle-diff.sh runs it.
set -euo pipefail
export LC_ALL=C TZ=UTC
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "${here}/../../.." && pwd)
ORACLE_SHA=98aa9252ff77825a8a797d388b449931fa9a234e

put() { mkdir -p "$(dirname "$1")"; printf '%b' "$2" > "$1"; at "$NEW" "$1"; }
at() { touch -t "$1" "${@:2}"; }   # [[CC]YY]MMDDhhmm[.SS], in UTC
NEW=202103040506.07  # after the run of the three groups' folders

for g in a b c; do mkdir -p "parts/$g"; done
(cd parts/a && bash "${here}/summary-collect.sh")
(cd parts/b && bash "${here}/summary-clinical.sh")
(cd parts/c && bash "${here}/summary-collect-c.sh")

# merge DEST SRC...: the folders' files, the later folder's file winning;
# times and symlinks kept.
merge() {
  local dest=$1 src
  shift
  mkdir -p "$dest"
  for src in "$@"; do cp -pPR "parts/${src}/." "$dest"/; done
}

mkdir -p assets py pymods
for g in a b c; do cp -pR "parts/${g}/assets/." assets/; done
cp -pR parts/b/pymods/. pymods/
cp -pR parts/c/py/. py/
git -C "$root" cat-file -e "${ORACLE_SHA}^{commit}" 2> /dev/null \
  || git -C "$root" fetch -q --depth=1 origin "$ORACLE_SHA"
git -C "$root" show "${ORACLE_SHA}:bin/run_delta.py" > py/run_delta.py

# all: every section of the three groups, the bash names.
merge all a/full b/cb-full c/full
# nf: the Nextflow names, steps failed and from an earlier run.
merge nf a/nf b/cb-nf c/nf
# bad: an unreadable file in sections of each group.
merge bad a/bad-a b/cb-bad c/bad1
# failed: steps that failed and wrote nothing, under not_run.
merge failed a/failed b/cb-failed
mkdir -p empty

# history: all, with three earlier runs kept. The newest has no summary
# (a half-copied folder) and is passed over; the next is compared: its
# ClinVar, CPIC, pypgx and PRS results and its section states differ, and
# so do the inputs of its manifest. A database re-check newer than the
# run's manifest.
merge history a/full b/cb-full c/full
put history/run_manifest.tsv "run\twritten_utc\t2021-03-04 05:06:07 UTC\nimage\tDEEPVARIANT_IMAGE\tgoogle/deepvariant:1.10.0\ndata\tclinvar_file_date\t2021-03-01\n"
put history/history/20200101T000000/summary.json '{"schema_version": 2, "sample": "S", "generated_utc": "2019-12-01 00:00 UTC", "sections": {}}'
put history/history/20210101T000000/summary.json '{"schema_version": 2, "sample": "S", "generated_utc": "2021-01-01 00:00 UTC",
 "sections": {
  "coverage": {"title": "Coverage (mosdepth)", "state": "missing", "data": {}},
  "clinvar": {"title": "ClinVar screen", "state": "ok", "data": {"hits": [
    {"chrom": "chr1", "pos": "45331556", "ref": "C", "alt": "T", "gene": "MUTYH", "significance": "Likely_pathogenic", "stars": 1, "review_status": "x", "genotype": "1/1"},
    {"chrom": "chr9", "pos": "1", "ref": "A", "alt": "G", "gene": "X", "significance": "Pathogenic", "stars": 2}]}},
  "cpic": {"title": "CPIC", "state": "ok", "data": {"genes": [{"gene": "CYP2D6", "diplotype": "*1/*1", "phenotype": "Normal"}]}},
  "pypgx": {"title": "pypgx", "state": "ok", "data": {"cyp2d6": "*1/*1"}},
  "prs": {"title": "PRS", "state": "ok", "data": {"scores": [{"pgs_id": "PGS000001", "condition": "x", "percentile": 10}]}},
  "mei": {"title": "Mobile element insertions (MEGAnE)", "state": "ok", "data": {}}},
 "cyp2d6": {"consensus": {"result": "*1/*1"}}}\n'
put history/history/20210101T000000/run_manifest.tsv "run\twritten_utc\t2021-01-01 00:00:00 UTC\nimage\tDEEPVARIANT_IMAGE\tgoogle/deepvariant:1.9.0\ndata\tclinvar_file_date\t2020-12-01\n"
mkdir -p history/history/20220101T000000
put history/history/20220101T000000/run_manifest.tsv "run\twritten_utc\t2022-01-01 00:00:00 UTC\n"
put history/recheck/state.json '{"checked_utc": "2021-03-05 00:00:00 UTC", "items": [
 {"item": "clinvar", "label": "ClinVar", "class": "newer-release", "ran_with": "2021-03-01", "upstream": "2021-03-04", "action": "rerun step 06", "extra": 1},
 {"item": "vep", "label": "VEP", "class": "no-newer-release"}, 7],
 "affected_steps": ["06"], "rerun": {"command": "run-all.sh"}}\n'

# history-bad: the only kept run has a summary that is no JSON, so nothing
# is compared; a re-check older than the run's manifest is of an earlier
# run and left out.
mkdir -p history-bad
put history-bad/run_manifest.tsv "run\twritten_utc\t2021-03-04 05:06:07 UTC\n"
put history-bad/history/20210101T000000/summary.json '{"schema_version": 2,'
put history-bad/recheck/state.json '{"checked_utc": "2021-03-01 00:00:00 UTC", "items": []}\n'

# inf: a telomere table whose read length is infinite, beside the
# coverage, step 33 and repeat sections. The Python collector stopped with
# an uncaught OverflowError and wrote no summary; og writes it with the
# telomere section unreadable (a golden).
mkdir -p inf
for f in mosdepth/S.mosdepth.summary.txt qc/S_sample_qc.tsv logs/run_status.tsv; do
  mkdir -p "inf/$(dirname "$f")" && cp -p "parts/a/full/$f" "inf/$f"
done
cp -pR parts/c/full/expansion_hunter inf/
put inf/telomere/S/S/S_summary.tsv "PID\tsample\ttotal_reads\tread_lengths\tintratel_reads\ttel_content\nS\ttumor\t1000\tinf\t40\t1.5\n"

rm -rf parts
