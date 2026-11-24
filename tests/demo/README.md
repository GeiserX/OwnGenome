# tests/demo

Tools for the docs pictures and the live-site demo. Every person here is invented: no real sample is read.

- `og dev demo-sample --genome-dir DIR` writes DIR/DEMO-001, the invented sample of the pictures, with every file steps 24, 27 and 28 read. Each table a step writes with og is written by that og verb, fed made-up inputs. The code is in `rust/crates/og-dev/src/demo/`, one file per report chapter.
  - `--sample DEMO-002 --sex female` and `--sample DEMO-003` write DEMO-001's invented partner and son for the live-site family page. DEMO-002 carries an FMR1 premutation, a CHIP call and a mosaic loss. DEMO-003 carries a class 4 deletion and two planted results to act on: a pathogenic ClinVar hit in KCNH2 at chr7:150960000 and a class 5 deletion over it.
  - `--seed N` changes the random values; the same seed writes the same bytes.
  - `og dev demo-sample check --genome-dir DIR` reads DIR/DEMO-001/summary.json and stops when a section is not ok, so a new report card that the demo does not feed yet is caught before a picture says "Not run".
- `retake-screenshots.sh WORK [PICTURES]` writes DEMO-001 with `bin/og`, runs steps 36, 27, 24 and 28 on it, runs the check, then draws and checks the pictures with `og dev pictures take` and `og dev pictures check` (in `rust/crates/og-dev/src/pictures.rs`), and checks with `og dev pictures same-sample` that the live-site demo shows the same Total variants. The demo-pictures workflow runs it. See docs/testing.md, "Retaking the report pictures".
- `og live demo --out DIR` builds the public site with three invented people, their reports, and a replay of DEMO-001. Run `og dev demo-events` to regenerate the replay log when the step catalogue changes.
