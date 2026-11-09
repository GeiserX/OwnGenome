# Direction A: Journal

Each person's genome reads like a long feature in a science journal. A serif headline and a lead paragraph say what matters in plain words. Numbered sections follow, with figures inside the text and caveats and evidence in the margin. A family is a report in chapters. A run in progress is the same report, with each section filling in as its steps finish.

## Navigation model

```
People and runs (home)
├── Being written: HG002  ──►  live run = the report, sections filling in place
├── Ruiz family  ──►  family report, chapters I–V + appendix
│     ├── Elena  ──►  person report, sections 1–10  ──►  finding page (e.g. BRCA2)
│     ├── Marc   ──►  person report (phone: one column, tab bar)
│     ├── Lía    ──►  person report (child: adult-onset results closed)
│     └── Carmen    (in the tree, no genome: dashed, "add a genome")
└── Every run (log for whoever runs the server)
```

- **One spine, three depths.** Family, then person, then finding, and the breadcrumb always shows where you are. You can go sideways at every depth: the person page's table of contents lists the other family members, a finding links to related findings (BRCA2 links to tamoxifen and CYP2D6, to the breast score, and to the family), and carrier results link to the couple screen.
- **Every page uses the same numbers.** The sections are 1 To act on, 2 Medicines, 3 For your family, 4 Worth a look, 5 Can I trust this?, 6 Origins, 7 Traits and scores, 8 Assessed, no finding, 9 What it cannot see, 10 Sources. The live run uses the same numbers, so a section being written sits where it will sit in the finished report.
- **Experts never leave the page.** The evidence sits in the margin of each finding (class, stars, gnomAD, reads, call quality). The finding page opens the full layer: ACMG codes, the frequency strip, a sketch of the reads with a link to the genome viewer, callers, and sources.

## What each page answers first

| Screen | First question it answers |
|---|---|
| 01 Home | "Does anyone have something to act on, and what is running?" The headline says it in a sentence, a "needs attention" list covers everyone, and the running run comes first. |
| 02 Live run | "Is it moving, when will it end, is anything wrong?" 31 of 54, about 19 h left, "nothing needs you". Then what has been found so far, written into its section. |
| 03 Person (Elena) | "Is there anything serious?" The lead paragraph names the condition (breast and ovarian cancer, BRCA2), then a row of counts by level and the data-quality line. |
| 04 Finding (BRCA2) | "What do I do about it?" A numbered plan sits just under the explanation, with the evidence below it. |
| 05 Family (Ruiz) | "What does it mean for our children?" The lead and chapter II: 1 in 4 for cystic fibrosis, drawn as four equally likely children, with Lía placed among them. |
| 06 Mobile (Marc) | "Which medicines should I watch out for?" After the lead, the medicine list comes first: drug, a verdict chip (Avoid, 20 mg max) and the gene as the reason. It is ready to show a pharmacist. |

## Five key design decisions

1. **The lead paragraph is the summary.** Every report opens with three or four plain sentences that name the conditions and drugs, not gene counts. A reader who stops there still knows what matters. The row of counts by level below it is the index.
2. **The margin holds caveats and evidence.** Each finding's caveat (first sentence visible, "More" for the rest), its evidence numbers and its definitions sit in a 260 px margin column beside it. Caveats stay visible without breaking the reading line, and experts get the numbers without a second page.
3. **Figures sit in the text, and each one does one job.** The figures are a chromosome map of the findings, the FMR1 repeat with its bands, small multiples for every repeat locus, PRS dots on a 0–100 strip with no coloured zones, a tall segmented ancestry bar with a chromosome painting, the BRCA2 protein with the cut marked, the frequency on a log strip, a relatedness dumbbell, the four-children grid, the parent-coloured chromosome map, and a step timeline. Charts use the sans, headings the serif. Status colours always come with a glyph and a label: diamond, pill, triangle, two rings, open circle, check.
4. **Silence is never an all-clear.** Every page has a section 8, "Assessed, no finding", that names each step, what it looked at and what it does not see, plus section 9 "What it cannot see". On Marc's phone, section 1 reads "No finding in what these steps assessed" and lists what those steps miss. Steps that failed or were skipped say so where it matters ("MUC1 kidney gene: not assessed, depth 41×"). APOE and an adult-onset result about a child stay closed behind an explicit "Show".
5. **Colour follows the person, only where people are compared.** Orange is Elena, blue Marc, aqua Lía, on the family page, the home tiles and the run log. Elena's own report draws no person colours, so the ancestry figures can use the categorical slots: no screen shows both sets. Green is used only for data and run states, never for a health result. Brass (the editorial accent) marks section numbers, links and "being written", and never carries a level.

## States covered

- **In progress:** sections stamped *Written 10:14*, *Being written* or *Waiting*, with grey proof lines and the step each one waits for. New findings carry a brass "New at 12:58" band.
- **Empty:** "No finding in what these steps assessed" plus what those steps do not see.
- **Failed or skipped:** an amber "Not assessed" card with the reason and a rerun action. Skipped steps are listed with plain reasons ("needs PacBio HiFi", "no symptoms given").
- **Sensitive:** a veiled card with an explicit button for APOE and for Lía's BRCA2 result.
- **Out of date:** "ClinVar 2026-11-06 is out" shows once in the sources and on the cards it would change (carrier list, couple screen). It is not shown on BRCA2, which it does not change.

## Files

- `screens/NN-screen.png`: desktop 1440 wide, full page. The mobile shot is 390 wide at 2x, at its full scroll height so the tab bar sits at the bottom.
- `html/NN-screen.html`: self-contained (CSS and JS inlined, system fonts only, no network). Built from `src/` by `node build.js`.
- `today/`: the live site before this design, for comparison.
- `BRIEF.md`: the data inventory and page list every design direction worked from; `dataset.json`: the invented family and findings the screens show.

## 07 Alt person (Elena, cover variant)

A bolder take on the same report, for comparison. The cover is a magazine spread: the name at 132 px, a one-sentence deck that names the condition, and a table of contents where each line already says what that section found. Below it, Elena's genome is drawn end to end as one band, with each finding pinned and labelled by condition. A five-column index lists the findings by name under each level, not as counts. The one finding to act on gets a feature spread with a pull quote, the protein cut and a four-step plan. Everything else becomes a plate of six figures (A to F), and section 8 is a ledger where each "no finding" line names what it cannot see.
