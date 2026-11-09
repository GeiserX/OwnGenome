# OwnGenome live site: design brief

Every designer works from this file and from `dataset.json` beside it. Four directions will show the same people and the same findings, so they can be compared on design alone.

The owner's verdict on today's site is that the UI and the UX are both bad. You should be able to see each family run and each person, with a report per person that you read on the page (not a PDF), and every result we found should be easy to see. He wants each page to look inventive and professional.

## 0. Hard rules

- **Dark theme only.** Screenshots use `colorScheme: 'dark'`.
- **Synthetic data only.** Use `dataset.json`: the Ruiz family (Elena, Marc, their daughter Lía, and Carmen, who is in the pedigree but not sequenced) plus a benchmark run of the public GIAB sample HG002. Never use or look for a real person's genome.
- **Static, self-contained HTML.** HTML, CSS and a little vanilla JS per screen. No CDN, no web-font downloads, no build step.
- **System fonts only.**
  - Display and headings: `ui-serif, "New York", "Iowan Old Style", Palatino, Georgia, serif`. The wordmark stays serif.
  - UI text: `system-ui, -apple-system, "Segoe UI", sans-serif`, with `font-variant-numeric: tabular-nums` on table columns and axis ticks.
  - Variant ids, HGVS and coordinates: `ui-monospace, "SF Mono", Menlo, monospace`.
  - Charts, numbers and figures stay in the sans. (The dataviz guidance says "no serif anywhere". The owner wants a serif for headings, so the serif goes in page headings only, never inside a chart.)
- **Screenshots.** Playwright, headless Chromium. Desktop is 1440x900 with `fullPage: true`. Mobile is 390x844 at `deviceScaleFactor: 2`. Name files `NN-screen.png`.
- **The site's own honesty rules stay** (`docs/live-site.md`; they are tested in `tests/test_live_results.py`):
  - **Never an all-clear.** No page says "normal", "fine", "healthy" or "nothing found" about a person's health. An empty section reads *"No finding in what this step assessed"*, followed by its caveat and what it does not see.
  - **APOE stays behind a gate** that the person opens on purpose.
  - **Polygenic scores show only with an ancestry-matched percentile.**
  - **Every card carries its caveat.** Its first sentence is visible, and the rest sits behind a toggle.
  - **A newer database release** shows on the cards it would change.
  - **A step that failed, was skipped or has not run** says so. It is never left blank.

## 1. Data inventory

The site has to show the result families below. The report already groups them into seven chapters, each named by a question (`bin/render_report.py` CHAPTERS). The tables keep that grouping. In each table:

- **One result** lists the fields a single result carries.
- **Read it as** is the plain-language way to read it.
- **Levels** uses the shared scale defined in section 4.3: *act, medicine, review, carrier, info, trust*.

### 1.1 Quality: "Can I trust this data?"

| Result | One result | Read it as | Levels | Caveats |
|---|---|---|---|---|
| Coverage (16b mosdepth) | mean depth, % of genome at 10x or more | about 30x is what was paid for; below 15x, small-variant calls lose accuracy | trust: ok, or warn below 15x | none |
| Sex check (16 indexcov, 33 somalier) | declared sex, sex from reads, chrX/chrY copy numbers | a mismatch means a swapped or mislabelled sample, or a real sex-chromosome difference; read nothing else until it is resolved | trust: ok / **stop** | none |
| Contamination (33 VerifyBamID2) | FREEMIX | above 0.03, someone else's DNA is mixed in | trust: ok / warn | none |
| Variant calling (03) | PASS variants (4.5 to 5.5 M), Ti/Tv (2.0 to 2.1) | a sanity check | trust | none |
| Callability (16b) | genes below 90% at 10x MAPQ 20, each with its %; genes "not on this reference" | **a "no finding" in these genes says nothing** | trust: warn list | SMN1, CYP21A2, PMS2, HBA, GSTT1 and others |
| Run health (live) | RAM, disk, retries, memory kills, database re-check | this describes the machine, not the person | trust | none |
| MultiQC (28) | a link out | expert only | none | none |

### 1.2 Clinical: "Is there anything to take to a doctor?"

| Result | One result | Read it as | Levels | Caveats |
|---|---|---|---|---|
| ClinVar screen (06) | gene, chrom:pos ref>alt, genotype (het/hom), CLNSIG (P/LP), review stars 0-4, gnomAD AF, ClinVar id, condition | 0-10 hits is typical. Most are recessive carriers. **Dominant gene, or two hits in a recessive gene = attention.** 2+ stars are more reliable | act / carrier / info | A frequency above 1% in gnomAD is unlikely to be high-penetrance. An unknown frequency has not been checked. Newer ClinVar may change the list |
| Secondary findings (23, ACMG SF v3.3, 84 genes) | gene, variant, the gene's reporting rule (dominant, two-hit recessive, HFE homozygous C282Y only, X-linked), TTN and LMNA notes | an unrequested but possibly serious, actionable finding (heart, cancer, cholesterol, metabolic) | act | one heterozygous MUTYH is listed apart as not reported |
| Cancer predisposition (17 CPSR) | per variant: tier/ACMG class (P, LP, VUS, LB, B), gene, consequence, evidence codes; counts per class | only P/LP need attention; **VUS = treat as not tested** | act (P/LP); VUS shown as counts only | the CPSR bundle's ClinVar release |
| Clinical filter (23) | about 300 rows: gene, impact (HIGH/MOD/LOW), CADD, REVEL, AlphaMissense, SpliceAI, LOEUF | expert worklist; read constrained genes and homozygous rows first | info (expert) | "a list to review with a clinician" |
| Variant prioritisation (31 slivar) | 4 groups; compound-het pairs with phase (trans/cis/unknown) | expert; most pairs are noise | info | phase is mostly unknown |
| Symptom ranking (48 Exomiser, opt-in) | top 10: rank, gene, score, mode of inheritance, variants, ACMG class; the HPO terms used | a ranking to take to a geneticist | review | not a diagnosis; an empty list rules nothing out |
| Structural variants (04/15/18/19/22/05) | type (DEL, DUP, INV, BND, INS), position, length, genotype, supporting callers, duphold DHFFC/DHBFC, AnnotSV class 1-5, genes, CDS overlap %, gnomAD-SV AF; counts per caller | 5-10k SVs is normal. **Only class 4-5 in a disease gene matters.** Single-caller BND calls are mostly artefacts | act (class 4-5) / info | short reads miss insertions and inversions |
| Mobile elements (50, opt-in) | 3 catalogued sites (FKTN, TAF1, BBS1): seen / not seen / not callable; about 1,200 Alu/L1/SVA calls | "not seen" is not proof of absence | carrier / review | none |
| Repeat expansions (09/09b) | locus, gene, two allele counts with 95% CI, Stranger status (normal, pre_mutation, full_mutation), normal max, pathogenic min, sized or flanking-only | HTT, FMR1, C9orf72, DMPK and ATXN1 first. **Intermediate = no disease, may expand in children. Full mutation = see a clinical geneticist.** An unsized locus is never shown as normal | act (full) / review (pre or intermediate) / info | short-read sizing fails on long expansions; C9orf72 thresholds vary between labs |
| MUC1 VNTR (52, opt-in) | positive / possible / negative / not assessed, sensitivity at this depth | ADTKD-MUC1, a kidney disease | act (positive) / review (possible) | needs confirmation (SNaPshot) |
| Mito variants (20) | site, allele fraction (heteroplasmy), ClinVar class, gnomAD AF; poly-C length variants listed apart | above 95% = haplogroup marker; 10-95% = look it up; m.3243A>G or m.8344A>G = specialist | act / review / info | blood under-represents heteroplasmy in other tissues |
| Lp(a) (46, opt-in, plus 42) | KIV-2 copies with 95% CI, rs10455872, rs3798220, flag "genetically high Lp(a)" | fewer copies mean higher Lp(a); **get one blood Lp(a) measurement** | review | a research estimate |
| CHIP (49, opt-in) | gene, protein change, VAF, alt reads/depth, class; per-gene detection-limit VAF | an acquired blood clone; common with age; most people with it never get a blood cancer | review / info | at 30x only VAF ≥ about 10-14% is visible, so a clean result does not rule CHIP out |
| Mosaic loss of Y/X (33) | % of cells with 95% range | becomes common with age | review / info | none |
| Mosaic chromosomal alterations (51 MoChA, opt-in) | chrom, start-end, Mb, arm, type (loss, gain, CN-LOH), cell fraction, reported y/n with reason | large changes carried by a share of the blood cells | review | detection floor in cell fraction |
| Somatic (29, experimental) | tumour-only Mutect2 calls | expert only; mostly germline false positives | info (expert) | highest false-positive rate in the pipeline |
| Viral (47, opt-in) | per virus: reads, % of the genome covered, depth relative to human; **integrated HHV-6** yes/no | iciHHV-6 makes a blood HHV-6 PCR read high. **Not a test for infection** | review (iciHHV-6) / info | reagent contamination |

### 1.3 Carrier: "What could I pass on to a child?"

| Result | One result | Read it as | Levels | Caveats |
|---|---|---|---|---|
| Carrier list (06 x 3,292 genes) | gene, variant, het, stars, condition, inheritance (AR/XL); hypomorph and pseudodeficiency alleles set aside | **carrier ≠ disease**; matters when the partner carries the same gene | carrier | "two hits across callers": a deletion plus a variant gives *biallelic* or *possibly biallelic* (then act) |
| Couple screen (39) | per gene: result (shared / x_linked / one_partner), offspring risk (1 in 4, 1 in 2), each partner's finding, condition, lists (ACMG/ACOG/Orphanet), explanation; read-depth alleles (GJB6, CFTRdele2,3, HBB, AZF); SMN1 copies; "does not see" list | **shared = 1 in 4 per pregnancy** | carrier (shared is the headline) | a screen, not a diagnosis; a list of what it misses |
| SMN1/SMN2, paralogs (35, opt-in) | agCN, psCN, quality, PASS; HBA copy count; CYP21A2 deletion/conversion | 1 SMN1 copy = SMA carrier; 2 copies does not rule it out | carrier | trust only quality ≥ 20 + PASS |

### 1.4 Medicines: "How might I respond to medicines?"

| Result | One result | Read it as | Levels | Caveats |
|---|---|---|---|---|
| PharmCAT (07) + CPIC (27) | gene, diplotype, phenotype (Normal, Intermediate, Poor, Rapid, Ultrarapid; Decreased function), status (normal, non-normal, ambiguous, not called, unclassified), phenoconversion; per drug: CPIC recommendation, strength, guideline | **one row per drug you might take** is what a person needs; share with the prescriber | medicine | "not called" ≠ normal; confirm before a prescription changes |
| CYP2D6 consensus (21/32/36) | pypgx vs Cyrius diplotypes, depth check, what was passed on (or `indeterminate`) | drug guidance only when two callers agree | medicine | the hard gene |
| HLA (08) | 2 alleles per gene; drug alleles present/absent (B\*57:01 abacavir, B\*58:01 allopurinol, B\*15:02/A\*31:01 carbamazepine); disease alleles (B\*27, C\*06:02); coeliac DQ2.5/2.2/8 with cis/trans | information about alleles; carrying one does not mean the reaction or the disease will happen | medicine (drug allele present) / review / info | approximate typing; not for transplant |
| Locityper (53, opt-in) | second opinion at HLA-DRB1, DQB1, KIR, CYP2D6, FCGR | shown **beside** the owner step's call, never instead of it | expert | none |

### 1.5 Ancestry: "Where do my ancestors come from?"

| Result | One result | Read it as | Caveats |
|---|---|---|---|
| Ancestry projection (26) | nearest panel population, PCA position | similarity, not a label | none |
| Fine ancestry (26b, opt-in) | continental (5) and regional (17) fractions with bootstrap low-high; chromosome painting (2 haplotypes x 22 chromosomes) | fractions with error bars | Oceania is not in the panel |
| mt haplogroup (12) | haplogroup, quality, contamination check (haplocheck) | the maternal line | population associations only |
| Y haplogroup (37, 45) | haplogroup, YFull terminal branch, path, private SNPs, Y-STR | the paternal line (males only) | none |
| Archaic (43, opt-in) | Neanderthal % and Denisovan %, segment count, 10 longest with genes | about 2% Neanderthal is typical outside Africa | none |
| Ancient (44, opt-in) | 20 nearest ancient groups with date, place and n; PCA with "you" marked | closeness is not descent | groups with n of 1-2 are noisy |
| ROH (11) | segments ≥ 5 Mb, total Mb, per chromosome; UPD flag | ~45 Mb second cousins, ~180 Mb first cousins | centromere artefacts |

All ancestry results are level **info**.

### 1.6 Traits and scores: "What do studies say about my traits and risk scores?"

| Result | One result | Read it as | Caveats |
|---|---|---|---|
| Traits (42) | trait, gene, rsid, genotype, plain reading, PMIDs; rules (APOE, TAS2R38, APOL1, SERPINA1 Pi, HFE, GC, 9p21, 17q21.31) | fun plus a few health-adjacent rows (APOL1, Pi\*Z, HFE) | APOE behind a gate |
| PRS (25) | condition, PGS id, percentile, ancestry group, matched %, applies-to sex | **a percentile is not a probability**; scores are not comparable with each other | hidden without the ancestry panel ("raw score only") |
| Blood groups (41, opt-in) | per system: phenotype, genotype, confidence, note (rare high-prevalence antigen missing goes first) | a prediction, not a serological typing | none |
| Telomere (10) | tel_content (relative) | compare only within one platform; **not biological age** | no normal range: show no healthy band |

All of these are level **info**. One exception: an APOL1 high-risk genotype, homozygous Pi\*ZZ or homozygous HFE C282Y is **review**.

### 1.7 Family (step 38, opt-in), on the family page

| Result | One result | Read it as |
|---|---|---|
| Relationship check | pair, declared, expected, measured relatedness, opposite homozygotes %, sites, verdict | **trust:** a disagreement stops the family analysis (sample swap) |
| Shared segments (hap-ibd) | pair, stretches, shared cM, IBD2 cM, longest, fraction, degree | confirms degree; finds distant cousins |
| Chromosome map | per child: segments from each parent (and grandparent when sequenced), crossovers, discordant sites | "which half of the genome came from whom" |
| Inheritance of findings | per finding in a parent: did the child inherit it (from `<child>_inheritance.tsv`) | **"what did my child inherit?"** |
| De novo | chrom/pos, child genotype, allele balance, GQ/DP of the three, pass/low, gene | most are genotype errors |
| Compound het | gene, the 2 variants, in trans / in cis / unknown | in trans in a fitting gene = review |
| Couple screen | see 1.3 | none |
| Plan | analysis, members, runs/skipped/waiting, reason | turn it into plain language, never step numbers |

Other and exports: the raw-data export (40, 23andMe-format file), the HTML report, the CPSR/PharmCAT/MultiQC reports, and the genome viewer (igv.js) at any position ("See the reads").

## 2. What people want to know, ranked

### A person opening their own genome

1. **"Is there anything serious?"** Count and name the *act* findings first, in words, with the condition and not only the gene. Today the ClinVar card says "hits: 3" and shows gene, position and stars, but **never the condition name**.
2. **"What do I do about it?"** One sentence per finding (see a genetic counsellor; confirm with a clinical test; ask for a blood Lp(a) test; tell your prescriber), plus who else in the family it concerns.
3. **"Which medicines should I watch out for?"** A drug list (drug, what to do, strength), not a gene table. It should be printable and easy to show a pharmacist on a phone.
4. **"What could I pass on, and what did my child inherit?"** The carrier list, the couple's shared risk with the 1-in-4 maths, and per child which findings came from whom.
5. **"Can I trust this?"** One line: data quality ok, identity confirmed, no contamination, plus what this genome could not read.
6. **"What is just interesting?"** Ancestry, haplogroups, Neanderthal %, traits, blood group, scores.
7. **"Is the run done, and is it healthy?"** Progress, an ETA, what failed and whether that matters.
8. **"What did it not look at?"** The "leaves out" list, so silence is never read as an all-clear.

### An expert (the owner, a clinician, a bioinformatician)

1. Per finding: the evidence (ClinVar stars and review status, gnomAD AF, ACMG codes, predictors), the reads (depth, allele balance, GQ, a link to the viewer), and the callers that agree.
2. QC detail: coverage distribution, callability per gene, somalier/FREEMIX, Ti/Tv, the depth needed per opt-in step.
3. Provenance: database releases used and newer ones available, tool versions, step state, retries, files.
4. Raw tables: every row sortable and searchable (the clinical filter's ~300 rows, all SVs, all STR loci, the slivar pairs).
5. Family: relatedness numbers, phase basis, de novo metrics, segment tables.

The layout follows from this: **the layperson's answer is the page, and the expert's evidence opens from it** (a disclosure, a side panel or a detail page). Experts are never sent to another page.

## 3. Today's site, and its problems

Screenshots (dark, captured from the demo replay, which holds one person and a two-member family):

- `current-01-home.png`: the home page
- `current-02-person.png`: a person, "Look first" tab
- `current-03-person-doctor.png`: a person, "Doctor" tab (6,433 px tall)
- `current-04-person-progress.png`: Progress
- `current-05-family.png`: the family page
- `current-06-person-mobile.png`: a person on a phone

### Hierarchy

- Every card has the same weight. The Doctor tab stacks 14 identical cards over 6,400 px. A pathogenic ClinVar hit and "Mosaic chromosomal alterations: none" look alike, with only a thin amber left border to tell them apart.
- **"Look first" lists the viral card as "worth a look" because its headline says "no integrated HHV-6".** The flag fires on a negative result, so the one curated list mixes signal with noise.
- The page top spends its space on a 54-cell progress strip and the sample id, even after the run has finished. Nothing there says what was found.
- Headlines are counts, not meanings: "Pathogenic/Likely Pathogenic hits: 3." The conditions (what each gene does to a person) appear nowhere on the person page.

### Density and jargon

- Stat grids carry labels like "Low-confidence calls (GQ below 30 or DP below 10), not counted", "By review stars (4/3/2/1/0): 0 / 0 / 2 / 1 / 0", "P/LP", "agCN", "Class 4-5 SVs", "F_ROH: 0.0000".
- Empty values show as `.`, `unknown` or `not checked` in table cells (ClinVar "Call" column, "gnomAD AF unknown" on every row).
- Each card's caveat sits in small grey text with a "Full caveat" toggle. The same "newer ClinVar available" chip repeats on many cards.
- Progress's "Found so far" repeats "finished, see the report" about 40 times. The ROH line is a 60-word paragraph.

### Navigation

- The tab names (Trust, Doctor, Child, Medicines, Ancestors, Traits, Other) are terse and only partly match the report's own questions. Tabs without findings look identical to tabs with them; only three carry a badge.
- On mobile the tabs scroll off-screen ("Me…" is cut) and the first screen shows two links and nothing else.
- Home is a near-empty page: one family header, two small cards, about 80% black space. People are shown as sample ids (DEMO-001), not names, and no card summarises what was found.
- To get from the family page to "what did the child inherit", or from a carrier finding to the couple screen, you have to go through tabs by hand. Findings never cross-link.

### Empty and odd states

- The search box is disabled: "Search arrives with the report".
- The family page says "2 of 2 sequenced" while showing DEMO-002 as "No run on this server" (dashed), yet it prints DEMO-002's carrier data.
- Empty states are written in pipeline terms: "Analyses: none planned yet (step 38 writes family_plan.tsv).", "Not run yet (step 38).", "The pedigree names no child."
- The couple screen table says "none shared" under "Chance for a child" and `.` for a partner without a finding. Its "Why" column is a 50-word paragraph per row.
- The step map is clipped inside its scroll box. Run health is four similar warning boxes.

### Visual

- Generic dashboard cards on one blue-grey (#14181f/#1d232d) with a blue link colour, so there is no brand. The serif appears only in the wordmark and two headings.
- Data is almost never visualised. There is one chromosome dot chart, and STRs, PRS, PGx, ancestry and pedigree are all tables or text.

## 4. Proposal

### 4.1 Pages: what each must answer

1. **Home: people and runs** (`/`). *Who is here, what is going on, and does anyone have something to act on?*
   - Families as units: a pedigree thumbnail, then the members as person tiles with name, role, age, run state, and per-level counts (act, medicine, review, carrier).
   - Runs in progress show live progress and an ETA. The unsequenced relative (Carmen) is shown, but quietly.
   - The benchmark run (HG002) is listed apart.
   - A "needs attention" strip across all people: "Elena: 1 finding to take to a doctor."

2. **Live run** (`/runs/HG002-BENCH/` while running). *Is it moving, when will it end, is anything wrong, and what has it found so far?*
   - Steps done out of planned, elapsed time, ETA, what runs now (CPU and memory per task), a step graph or timeline, events, machine health, and "found so far" in plain words.
   - A finished run turns this page into the report and keeps Progress as a secondary view.

3. **Person report** (`/runs/RUIZ-ELENA/`). *The report as a page.*
   - Order: a summary (who, the trust line, a count by level), then **To act on**, **Medicines**, **For your family**, then **Explore** (ancestry, traits, scores, all sections).
   - Each finding is a readable block: what it is, what it means, what to do, who else it concerns, its caveat, and evidence on demand.
   - "No finding in what was assessed" sections collapse into one compact list that still names each section and its caveat.
   - Shows the leaves-out list, the database freshness, and the files and report links.

4. **Finding detail** (`/runs/RUIZ-ELENA/#f/F01`, a page or a side panel). *Everything about one result.*
   - The plain explanation, the evidence (classification, stars, review, gnomAD, predictors), reads (depth, allele balance, GQ, viewer link), inheritance in the family, related findings (BRCA2 links to tamoxifen/CYP2D6 and the PRS), the caveat in full, and sources.

5. **Family** (`/runs/RUIZ-FAMILY/`). *How are we related, what does each child carry from whom, and what does it mean for future children?*
   - A large pedigree with per-person finding markers.
   - The couple screen headline (shared CFTR, 1 in 4) with the risk drawn as a 2x2 grid of outcomes.
   - Lía's inheritance of findings: a table with parent colours.
   - The chromosome map (two rows per chromosome, coloured by parent).
   - Relationship checks (measured against expected), de novo, compound het, the plan in plain words, and what the screen does not see.

6. **Mobile person** (390 wide). *The same report, one column.*
   - A sticky summary, findings as stacked blocks, and the medicines list easy to show a pharmacist.
   - The section nav becomes a bottom sheet or a segmented control, not a sideways-scrolling tab row.

### 4.2 Shared dataset

`dataset.json` holds everything:

| Section | What it holds |
|---|---|
| `people` | 5 entries: Elena (38, F), Marc (41, M), Lía (6, F, a minor), Carmen (not sequenced), HG002 (public reference) |
| `runs` | three finished runs (Marc's has one failed opt-in step, 52 MUC1, at too low a depth), the family run with its plan, and HG002 running (31/54, 4 tasks live, events, found so far, machine health) |
| `families[0]` | relationships, shared segments, the couple screen (CFTR shared; PAH, HFE one partner; the FMR1 context row; read-depth alleles; SMN1; the does-not-see list), and the trio (inheritance of findings, de novo, compound het, chromosome-map counts and a chr7 example) |
| `findings` | 28 findings with `level`, `chapter`, `section`, plain-language `what_it_means`, `what_to_do`, `caveat` and the evidence fields |
| `person_panels` | STR loci, ROH, telomere, archaic, traits, blood, PRS lists, mosaic and viral, for the "Explore" parts |

The findings by level:

| Level | Count | Findings |
|---|---|---|
| act | 1 | Elena BRCA2 |
| medicine | 4 | Elena CYP2D6 IM; Marc CYP2C19 PM; Marc SLCO1B1; Lía CYP2C19 IM |
| review | 4 | Elena FMR1 intermediate; Lía FMR1 intermediate; Marc high Lp(a); Marc HLA-B\*27 |
| carrier | 6 | CFTR x3 (Elena, Marc, Lía, whose copy is inherited from Elena); PAH; HFE; the couple CFTR row |
| info | 10 | including APOE, which is gated |
| trust | 3 | none |

The stories the data is built to tell. Every direction must make each of them readable within a few seconds:

- **Elena has one finding to act on (BRCA2).** It links to her CYP2D6 result (tamoxifen) and to "Lía did not inherit it". That last line is an adult-onset finding about a minor, so put it behind an explicit reveal.
- **Marc has no act finding but two medicine findings** (clopidogrel: avoid; simvastatin: 20 mg or less) and a heart cluster: high Lp(a), CAD score at the 88th percentile, SLCO1B1. Show them together.
- **The couple shares CFTR (1 in 4).** Lía inherited Elena's copy only, so she is a carrier and not affected.
- **FMR1 48 in Elena became 49 in Lía:** intermediate in both.
- **Lía is O RhD-negative from two RhD-positive parents.** Show the inheritance.
- **HG002 is mid-run:** one memory retry, four tasks running, early findings arriving.

### 4.3 One severity scale for every result

| Level | Words on screen | Colour (status palette, dark) | Icon | Rule |
|---|---|---|---|---|
| act | "Take to a doctor" | critical `#d03b3b` | filled diamond or ! | always pair with the condition name and "what to do" |
| medicine | "Tell your prescriber" | serious `#ec835a` | pill or Rx | drug-first wording |
| review | "Worth a look" | warning `#fab219` (dark text on fill) | ▲ | a conversation or a measurement, not a diagnosis |
| carrier | "For family planning" | a neutral accent (e.g. violet `#9085e9` used as an outline, never a fill that reads as alarm) | two-circle glyph | always say "not affected" |
| info | "For information" | no colour: secondary ink | ○ | never "normal" or "fine" |
| trust | "Data quality" | good `#0ca30c` **only for data and run states**, warn/critical for problems | ✓ / ! | green is never used for a health result |

Status colours always come with an icon and a label, never colour alone. They are never reused as series colours.

### 4.4 Visual system starting points

Each direction may diverge from these, but stays dark and keeps these roles.

- **Surfaces.** Page `#0f1113` to `#121417`; raised `#17191d`; hairline `rgba(255,255,255,.08)`.
- **Ink.** Primary `#f0efec`; secondary `#c3c2b7`; muted `#898781`.
- **Person identity colours.** Validated with the dataviz validator, dark, all pairs, surface `#17191d`: all checks pass.
  - Marc: blue `#3987e5`
  - Elena: orange `#d95926`
  - Lía: aqua `#199e70`
  - Use them for the pedigree, the chromosome map, relatedness and family tables. **Colour follows the person** on every page.
  - On a person page, ancestry charts may reuse the categorical slots. Never draw person colours and ancestry colours on one screen.
- **Chart forms**, chosen by each result's job (dataviz form heuristic):
  - **Findings summary:** stat tiles or a level-count row, not a donut.
  - **PRS:** a dot on a 0-100 percentile strip, one row per condition, with the condition's ancestry group named. No red/green zones, because a percentile is not a risk class.
  - **STR:** per locus, the two allele dots on a count axis, with labelled bands (normal, intermediate, premutation, full) in muted tones plus text labels. An unsized allele is drawn as an open-ended bar from its lower bound.
  - **Lp(a):** the KIV-2 count with its CI on a strip, plus the two variants as chips. **CHIP:** a meter of the detection limit. **Telomere:** a plain number only, because it has no normal range.
  - **Ancestry:** a tall segmented bar (≤ 6 parts, with "Other"), plus a 22-chromosome painting with 2 haplotypes each.
  - **Variant locations:** a chromosome ideogram with a marker per finding.
  - **ROH:** a per-chromosome bar of Mb in runs ≥ 5 Mb.
  - **Mito:** heteroplasmy dots on a 0-100% axis.
  - **Couple risk:** a 2x2 grid of the four equally likely children (1 affected, 2 carriers, 1 neither), labelled.
  - **Relatedness:** a dumbbell of expected against measured per pair. **Chromosome map:** paternal row over maternal row, in person colours.
  - **PGx:** a table or cards, drug first, with gene → phenotype as the reason. Not a chart.
  - **Live run:** progress as a hero figure (31 of 54), a step timeline or dependency graph coloured by state (waiting, running, done, failed, skipped) with labels, and CPU/memory sparklines per running task.
- **Chart rules.** Thin marks, recessive grid, one axis, direct labels, and a hover tooltip per mark. Every chart has a table view or its numbers in text beside it. Run `scripts/validate_palette.js` on any new categorical set.

### 4.5 Copy rules

- Names, not sample ids ("Elena Ruiz"; the run id goes in small print).
- The condition name before the gene: "Hereditary breast and ovarian cancer (BRCA2)".
- No step numbers or file names in body text. They belong in the evidence and provenance layer.
- Define a term the first time it appears (an inline definition or a tooltip): carrier, heterozygous, percentile, VUS, heteroplasmy.
- Keep caveats short, specific and visible. Never put them in a grey footnote as a disclaimer.
- An empty cell reads "none found", "not assessed" or "not applicable". It is never `.`.

### 4.6 What each direction delivers

Six screens, using `dataset.json`, saved as `NN-screen.png`:

1. home (desktop)
2. live run, HG002 (desktop)
3. person report, Elena (desktop)
4. finding detail, Elena BRCA2 (desktop)
5. family, Ruiz (desktop)
6. person report, Marc (mobile)

Plus one line per screen saying which question in section 2 it answers first.
