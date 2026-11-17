# User journeys

This page records where a person gets stuck on each OwnGenome surface. Each row is one journey: a goal, a surface and the outcome of each walk. When a fix lands, walk its journey again and add the new outcome with its date and commit.

The first run took place on 2026-11-17, on commit `65095e97` (main). It used only invented demo data (DEMO-001, DEMO-002, DEMO-003 and their family) and public test data. No real genome was used.

## How the run worked

- Agents played each persona. They started from a public entry point (the README, the docs site or the machine-readable docs index) and never read the source.
- A new-user journey had 20 actions. An advanced, self-hosting or agent journey had 40.
- Each new-user walk ran twice, in separate fresh sandboxes, to separate a real stuck point from chance.
- A scorer read each transcript and gave every journey an outcome. A separate verifier then tried to reproduce each finding in its own fresh browser or shell.
- Each kit got a planted control, a blocker put there on purpose, to prove the walk can report a stuck point.
- Sweeps checked the docs site and the demo at 320, 390, 820 and 1440 px, with failed and slow data, keyboard use and colour contrast.

Of 138 candidate findings from the walks, the verifier reproduced 78. Three of them showed only some of the time. It found 25 caused by the test environment alone and could not reproduce 33. Two more were the planted controls. The sweeps gave 15 findings, and the verifier reproduced 14.

### Outcomes

| Outcome | Meaning |
| --- | --- |
| done | The person reached the goal without help. |
| done with a minor issue | The goal was reached, with friction a person would notice. |
| done with help | The goal was reached, but the walker did something a person might not. Here it was mostly a click on a link or control that was below the screen in its last screenshot. A person would scroll first. |
| failed | The goal was not reached. The blocker column says why. |
| not reached | The walk never got to this journey. |

Severity says what the person loses. A blocker stops the goal. A major finding gives a wrong or missing result. A minor finding costs time. A cosmetic finding costs nothing but polish.

## Surfaces

| Surface | Entry point | Walked |
| --- | --- | --- |
| Docs site | https://geiserx.github.io/OwnGenome/ | desktop and phone |
| README | https://github.com/GeiserX/OwnGenome | desktop, phone and shell |
| Live demo with invented runs and a family | https://geiserx.github.io/OwnGenome/demo/ | desktop and phone |
| HTML report and text exports, opened from the demo | the demo | desktop and phone |
| Setup and run scripts, and the `og` command | the README install block | shell |
| Local live site with a password | [Watch a Run Live](live-site.md) | not reached |
| Standalone read viewer | [Genome Viewer](viewer.md) | not reached |
| Nextflow route | [Nextflow Execution](nextflow.md) | not reached |
| Machine-readable docs index | https://geiserx.github.io/OwnGenome/llms.txt | agent |
| Public test data releases | [Testing](testing.md) | shell |

OwnGenome has no hosted account, email sign-up, mobile app or MCP server, so no journey covers them.

## Personas

| Id | Persona | Kit | Knows | Actions per journey |
| --- | --- | --- | --- | --- |
| N01 | New user choosing a route | desktop browser | the public introduction | 20 |
| N02 | New user choosing a route | phone browser | the public introduction | 20 |
| N03 | New user trying the demo | desktop browser | nothing yet | 20 |
| N04 | New user trying the demo | phone browser | nothing yet | 20 |
| N05 | New user installing from the README | shell | the README | 20 |
| N06 | New user following the quick test | docs and shell | the docs | 20 |
| A01 | Advanced user with public reads | shell | file formats and the docs | 40 |
| A02 | Advanced user with chip data | shell | the chip-data guide | 40 |
| A03 | Advanced user with clinical and population analyses | docs and shell | the analysis docs | 40 |
| A04 | Advanced user with optional analyses | docs and shell | the pipeline overview | 40 |
| H01 | Self-hoster preparing a local install | docs and shell | the published install docs | 40 |
| H02 | Self-hoster using the install on desktop | desktop browser | connection details and an invented password | 40 |
| H03 | Self-hoster using the install on a phone | phone browser | connection details and an invented password | 40 |
| H04 | Self-hoster keeping the install running | docs and shell | a copy of earlier outputs | 40 |
| G01 | Coding agent with the README and help | shell | the README, linked docs and help | 40 |
| G02 | Coding agent with the docs index | docs and shell | the machine-readable docs index | 40 |

The desktop browser was 1440 x 900. The phone browser emulated an iPhone 15 with touch. The shell kit was a clean home folder on an Apple Silicon Mac with 16 GB of memory.

## Journeys and outcomes

The first useful result is measured once per walk, from its first action. A walk covers several journeys, so its rows share the figure. Where only docs or help appeared, the cell says so. For the twice-walked journeys, the cell gives the first walk, then the repeat.

The blocker codes point to [what this run could not check](#what-this-run-could-not-check).

| Journey | Persona | Surface | Outcome | Blocker | Date | Commit | First useful result |
| --- | --- | --- | --- | --- | --- | --- | --- |
| D01 | N01 | Published docs, desktop | done with a minor issue, then done with a minor issue | none | 2026-11-17 | 65095e97 | 18 s, then 16 s |
| D02 | N01 | README and docs, desktop | done with help, then done with help | none | 2026-11-17 | 65095e97 | 18 s, then 16 s |
| P01 | N02 | Published docs, phone | done with help, then done with help | none | 2026-11-17 | 65095e97 | 20 s, then 11 s |
| P02 | N02 | README and docs, phone | done with help, then failed | action budget on the repeat | 2026-11-17 | 65095e97 | 20 s, then 11 s |
| W01 | N03 | Demo found through the docs, desktop | done, then done with a minor issue | none | 2026-11-17 | 65095e97 | 37 s, then 29 s |
| W02 | N03 | Demo while a run is going, desktop | done with help, then failed | a reload resets the replay moment | 2026-11-17 | 65095e97 | 37 s, then 29 s |
| W03 | N03 | Demo family, desktop | done with help, then done with help | none | 2026-11-17 | 65095e97 | 37 s, then 29 s |
| W04 | N03 | Demo gene and variant location, desktop | done with help, then done with help | none | 2026-11-17 | 65095e97 | 37 s, then 29 s |
| W05 | N03 | Demo report files, desktop | done with help, then done with help | none | 2026-11-17 | 65095e97 | 37 s, then 29 s |
| T01 | N04 | Demo found through the docs, phone | done with help, then done with a minor issue | none | 2026-11-17 | 65095e97 | 1 min 41 s, then 1 min 30 s |
| T02 | N04 | Demo while a run is going, phone | done with help, then done with help | none | 2026-11-17 | 65095e97 | 1 min 41 s, then 1 min 30 s |
| T03 | N04 | Demo family, phone | done with help, then done with help | none | 2026-11-17 | 65095e97 | 1 min 41 s, then 1 min 30 s |
| T04 | N04 | Demo gene and variant location, phone | done with help, then done with help | none | 2026-11-17 | 65095e97 | 1 min 41 s, then 1 min 30 s |
| T05 | N04 | Demo report files, phone | done with help, then failed | no save dialog in the test phone browser, export files named as text | 2026-11-17 | 65095e97 | 1 min 41 s, then 1 min 30 s |
| I01 | N05 | README install, shell | failed, then failed | M, Docker memory check | 2026-11-17 | 65095e97 | 1 min 56 s, then 3 min 33 s (checks only) |
| I02 | N05 | First run, shell | failed, then failed | J | 2026-11-17 | 65095e97 | 1 min 56 s, then 3 min 33 s (checks only) |
| I03 | N05 | Validation and recovery, shell | failed, then failed | M | 2026-11-17 | 65095e97 | 1 min 56 s, then 3 min 33 s (checks only) |
| Q01 | N06 | Quick-test setup | failed, then failed | M, W, action budget | 2026-11-17 | 65095e97 | 40 s, then 25 s (docs only) |
| Q02 | N06 | Quick-test VCF jobs | failed, then failed | Q01 not ready | 2026-11-17 | 65095e97 | 40 s, then 25 s (docs only) |
| Q03 | N06 | Quick-test aligned reads | failed, then failed | Q01 not ready, W | 2026-11-17 | 65095e97 | 40 s, then 25 s (docs only) |
| R01 | A01 | Public test data, shell | failed | M | 2026-11-17 | 65095e97 | 4 min 1 s (download only) |
| R02 | A01 | Raw reads to a report, shell | failed | M, J | 2026-11-17 | 65095e97 | 4 min 1 s (download only) |
| R03 | A01 | BAM and CRAM input, shell | failed | M, J | 2026-11-17 | 65095e97 | 4 min 1 s (download only) |
| R04 | A01 | VCF and gVCF input, shell | failed | M, J | 2026-11-17 | 65095e97 | 4 min 1 s (download only) |
| R05 | A01 | Settings, restart and run state, shell | failed | M, J | 2026-11-17 | 65095e97 | 4 min 1 s (download only) |
| R06 | A01 | Lookup and archive, shell | failed | M (`og` never built) | 2026-11-17 | 65095e97 | 4 min 1 s (download only) |
| C01 | A02 | 23andMe chip import, shell | failed | R | 2026-11-17 | 65095e97 | none |
| C02 | A02 | MyHeritage chip import, shell | failed | R | 2026-11-17 | 65095e97 | none |
| C03 | A02 | AncestryDNA chip import, shell | not reached | R | 2026-11-17 | 65095e97 | none |
| C04 | A02 | Chip analysis and limits, shell | not reached | R | 2026-11-17 | 65095e97 | none |
| C05 | A02 | Raw-data export and imputation, shell | not reached | R | 2026-11-17 | 65095e97 | none |
| B01 | A03 | Clinical annotation | failed | R, M | 2026-11-17 | 65095e97 | 57 s (docs only) |
| B02 | A03 | Cancer and symptom ranking | failed | R, M | 2026-11-17 | 65095e97 | 57 s (docs only) |
| B03 | A03 | Medicines and HLA | failed | R, M | 2026-11-17 | 65095e97 | 57 s (docs only) |
| B04 | A03 | Structural variants and repeats | failed | R, M, J | 2026-11-17 | 65095e97 | 57 s (docs only) |
| B05 | A03 | Polygenic scores and ancestry | failed | R, M, J | 2026-11-17 | 65095e97 | 57 s (docs only) |
| B06 | A03 | Fine, archaic and ancient ancestry | failed | R, M, J | 2026-11-17 | 65095e97 | 57 s (docs only) |
| B07 | A03 | Mitochondria, lineage and telomeres | failed | R, M, J | 2026-11-17 | 65095e97 | 57 s (docs only) |
| B08 | A03 | Combined report and rerun changes | failed | R, M, J | 2026-11-17 | 65095e97 | 57 s (docs only) |
| X01 | A04 | Traits, blood groups and copy number | failed | R | 2026-11-17 | 65095e97 | none |
| X02 | A04 | Somatic and viral analyses | not reached | R | 2026-11-17 | 65095e97 | none |
| X03 | A04 | Difficult loci and mobile elements | not reached | R | 2026-11-17 | 65095e97 | none |
| X04 | A04 | Other callers and benchmarking | not reached | R | 2026-11-17 | 65095e97 | none |
| X05 | A04 | Family readiness | not reached | R | 2026-11-17 | 65095e97 | none |
| H01 | H01 | Local install with no data yet | failed | M | 2026-11-17 | 65095e97 | 42 s (docs), 4 min to a started site |
| H02 | H01 | Local example with results and reads | failed | M (`og` never built) | 2026-11-17 | 65095e97 | 42 s (docs), 4 min to a started site |
| H03 | H01 | Password and session length | failed | H02 not ready | 2026-11-17 | 65095e97 | 42 s (docs), 4 min to a started site |
| H05 | H01 | Standalone read viewer | failed | H02 not ready | 2026-11-17 | 65095e97 | 42 s (docs), 4 min to a started site |
| H04 | H01 | Live progress from a small run | failed | J | 2026-11-17 | 65095e97 | 42 s (docs), 4 min to a started site |
| U01 | H02 | Sign in, sign out and return, desktop | failed | H | 2026-11-17 | 65095e97 | 54 s (docs only) |
| U02 | H02 | Empty and partial results, desktop | failed | H | 2026-11-17 | 65095e97 | 54 s (docs only) |
| U03 | H02 | Live and standalone reports, desktop | failed | H | 2026-11-17 | 65095e97 | 54 s (docs only) |
| U04 | H02 | Read evidence and family results, desktop | failed | H | 2026-11-17 | 65095e97 | 54 s (docs only) |
| U05 | H02 | Session expiry, desktop | failed | H | 2026-11-17 | 65095e97 | 54 s (docs only) |
| V01 | H03 | Sign in, sign out and return, phone | failed | H | 2026-11-17 | 65095e97 | 55 s (docs only) |
| V02 | H03 | Empty and partial results, phone | failed | H | 2026-11-17 | 65095e97 | 55 s (docs only) |
| V03 | H03 | Live and standalone reports, phone | failed | H | 2026-11-17 | 65095e97 | 55 s (docs only) |
| V04 | H03 | Read evidence and family results, phone | failed | H | 2026-11-17 | 65095e97 | 55 s (docs only) |
| V05 | H03 | Session expiry, phone | failed | H | 2026-11-17 | 65095e97 | 55 s (docs only) |
| M01 | H04 | Settings and logs | failed | C | 2026-11-17 | 65095e97 | 54 s (docs only) |
| M02 | H04 | Database recheck and scheduling | failed | C | 2026-11-17 | 65095e97 | 54 s (docs only) |
| M03 | H04 | Backup and restore | failed | C, no documented procedure | 2026-11-17 | 65095e97 | 54 s (docs only) |
| M04 | H04 | Upgrade and proxy under a path | failed | C, no upgrade guide | 2026-11-17 | 65095e97 | 54 s (docs only) |
| M05 | H04 | Leaving and cleanup | failed | C | 2026-11-17 | 65095e97 | 54 s (docs only) |
| G01 | G01 | Agent install and first run | failed | M | 2026-11-17 | 65095e97 | 2 min 24 s (help only) |
| G02 | G01 | Agent public-data analysis | failed | M, G01 not ready | 2026-11-17 | 65095e97 | 2 min 24 s (help only) |
| G03 | G01 | Agent gene question and portable output | failed | G02 not ready | 2026-11-17 | 65095e97 | 2 min 24 s (help only) |
| G04 | G01 | Agent recovery and repeat | failed | M | 2026-11-17 | 65095e97 | 2 min 24 s (help only) |
| L01 | G02 | Agent picks a route from the docs index | done with help | none | 2026-11-17 | 65095e97 | 52 s |
| L02 | G02 | Agent runs the small VCF example | failed | Docker memory check, B | 2026-11-17 | 65095e97 | 52 s |
| L03 | G02 | Agent checks if results need updating | done with help | none | 2026-11-17 | 65095e97 | 52 s |
| F201 to F211 | various | Bad input, failed service, expired session, destructive recovery | not reached | round 2 not walked | not walked | not walked | not walked |
| F301 to F307 | various | Phone and agent parity, family results, remote server, real volume | not reached | round 3 not walked | not walked | not walked | not walked |

## Findings with a fix in progress

Each row is one root cause. Its fix is named by the pull request title. Most of these pull requests are not open yet.

| Fix | What a person hits | Severity | Status |
| --- | --- | --- | --- |
| fix(demo): replay the demo as one run so unstarted sections say waiting, not written | At the start of the demo replay, sections read Written while the page says no section is written yet. Midway, the medicines card says no step made it, beside Waiting for PharmCAT. The report gives a run date in October, the replay ends in January. | blocker | fix in progress |
| fix(live): stop the family page calling flagged cards 'worth a look' when the person page counts six results | The family page says 2 and 5 results are worth a look for the two parents. Their own pages and the overview say 6 and 9. | major | fix in progress |
| fix(live): keep a child's opened result on the family page instead of replacing it with the inheritance chapter | A child's last closed result shows for a moment, then the inheritance chapter replaces it. On a phone, closing it leaves the controls above the screen. A reloaded chapter link opens the family introduction. | major | fix in progress |
| fix(docs): make the Nextflow install and three other table commands copy and run as written | The Nextflow install command inside a docs table prints a backslash before the pipe, so the copied command does not run the installer. Three other table commands do the same. | major | https://github.com/GeiserX/OwnGenome/pull/493 |
| docs(quick-test): say the test VCF covers chr1-22 and weighs 120 MB, so the chr22 step makes sense | The quick test calls its download chr22 only, then full genome. It says about 250 MB, and the file is about 120 MB. | minor | https://github.com/GeiserX/OwnGenome/pull/491 |
| docs(live): link the demo, say how to stop the site, and say it needs its own host name behind a proxy | The demo address on the docs home page is plain text, not a link. The live-site guide does not say how to stop the site. It gives no way to serve it under a path behind a proxy. | major | https://github.com/GeiserX/OwnGenome/pull/492 |
| docs(install): give requirements first, with Homebrew before brew and the swap note on 16 GB | The Mac steps say `brew install bash` but never say to install Homebrew. The 16 GB minimum leaves out that a 16 GB machine needs swap for some tasks. | minor | fix in progress |
| fix(run-all): let the run be told which input to start from, and say when the gVCF is read | Two guides disagree on whether `run-all.sh` reads the gVCF that sits beside an existing VCF. | major | fix in progress |
| docs(ops): stop saying exit 137 always means out of memory, and that og needs a hand-built binary on a Mac | Two guides give different steps for running `og` on a Mac. | major | fix in progress |
| fix(live): label the reads drawing an illustration and say the demo has no reads before the click | See the reads opens a panel below the screen, and only there says the demo has no reads. | minor | fix in progress |
| feat(live): build the public demo with og, labelled as invented in every file it serves | The downloaded demo report and medicine text do not say the people are invented. Their version fields read unknown. A gene track says no gene beside a marked variant. The CPSR link in each demo report opens a missing page. | major | fix in progress |
| build(demo): take the demo pictures in Rust from the same sample the live demo serves | The report picture in the docs and the demo report give different variant totals for the same invented sample. | major | fix in progress |
| fix(report): show PGx incomplete calls, the CYP2D6 flags and contraindications in both reports | The CYP2D6 depth check says passed while every depth value under it reads not measured. | major | fix in progress |
| fix(pgx): credit PharmGKB in the drug guidance it supplies | One drug's guidance prints `&lt;` where it means the less-than sign. | cosmetic | fix in progress |
| fix(live): add skip links and finger-sized caveat, replay and timeline controls | On a phone, the step timeline scrolls inside the page and a tapped step's detail stays out of view. The More controls are 15 px high. | minor | fix in progress |
| fix(report): put the research-only line first, fit tables on a phone and add a skip link | Report cards and tables run past the screen at 320, 390 and 720 px, so a person scrolls sideways to read them. | major | fix in progress |
| fix(live): a person or home page that cannot load its data says so and offers a retry | When the demo data fails or is slow, a page goes blank or says nothing finished yet, with no retry. A failed adult page reads as a child's closed report. | minor | fix in progress |
| fix(live): stop the genome map's zoom buttons pushing the person page 20 px past a small phone's screen | Both adult person pages are 20 px wider than a 320 px screen. | minor | fix in progress |
| fix(live): stop offering 'Open in the report' for a section the report leaves out | Open in the report opens the report at the top, because the report has no such section. | minor | fix in progress |
| fix(report): make the Exomiser and trait annotation lines readable on the dark report | The small annotation lines in the report have a contrast of about 2.8 to 1. | minor | fix in progress |
| fix(docs): make the current page readable in the side menu and drop 122 empty tab stops from code blocks | In dark mode, the current page in the docs menu has a contrast of 3.1 to 1. Keyboard users pass 122 unnamed code-line links on one page. | minor | fix in progress |

## Findings that need a decision or another project

These need a product decision, a new procedure or a change in the docs theme. The last column says who should take each.

| What a person hits | What closes it | Who |
| --- | --- | --- |
| The docs ask for at least 16 GB of Docker memory, which is all of a 16 GB Mac. The setup check fails below 8 GB. A 16 GB Mac with Docker's default gets far too low. | One Docker memory figure in every doc and check, proven by a default run on a 16 GB Mac, or a higher stated minimum. | owner |
| There is no guide to upgrade an existing installation. The updating page is for maintainers. | An upgrade page, proven by upgrading from the previous release. | owner |
| There is no way documented to back up and restore results together with the live site's state. | A section that names every folder to keep, how to restore it and how to check the result, tested once. | owner |
| The cleanup steps prune every Docker image on the machine. Nothing says how to remove only OwnGenome. | An uninstall section for the pinned images, the data folders and the service, with the global prune as a last resort. | owner |
| The demo export card and the report name the export files and the sharing guide as plain text, not as links. | A decision on whether the site offers those files, and one rule for whether report notes link to docs pages. | owner |
| An exome has no single command. The exome guide lists the steps one by one. | An exome mode in `run-all.sh`, or a pointer to the Nextflow exome mode once it lands. | owner |
| A reload of the demo replay returns to the finished moment. | A decision on keeping the moment in the address. | owner |
| On a phone, docs tables and long commands scroll sideways, so the first column leaves the screen. | A phone rule for tables in the docs theme, checked by a walk at 393 px. | owner |
| On a page scroll, the docs site can show a blank frame or overlapping titles for under a second. | A recording in a real browser, then an issue for the docs theme or a fix in our stylesheet. | upstream (Zensical) |
| On a phone, the docs contents panel stays open over the section it jumped to. | A theme fix or option that closes the panel, checked at 393 px. | upstream (Zensical) |
| The code copy button on a phone is 22.5 px square. | A look at the theme's layout, or a theme change. | upstream (Zensical) |
| The footer credit is 4.32 to 1 in dark mode, just under 4.5 to 1. | One colour rule in our stylesheet, if wanted. | owner |
| With text-only enlargement, the demo finding title breaks into letter groups. Normal browser zoom works. | A decision on supporting text-only enlargement, then layouts sized in text units. | owner |

## What this run could not check

The test machines limited what the run could reach. These limits come from the environment, not from OwnGenome.

| Code | Limit | What would reach it |
| --- | --- | --- |
| J | No Java VM starts on the test Macs' macOS version (26.6.1). Every JDK tried crashed at start, so no Nextflow or pipeline run could finish. | Walk again on Linux, or on a Mac whose macOS version starts Java. |
| M | The container runtime there (Colima) shared only the home folder. The walk folders were on another disk, so most shell setup walks failed on bind mounts. This is an environment limit, not a product defect. | Put the walk folders under the home folder, or use Docker Desktop as the docs describe. |
| R | Three docs and shell walks had a docs reader with no browser, and it could not open the docs site. | Give those walkers a browser, as a person reads the docs in one. |
| H | The browser walks of the local site could not read the connection details. The local install walk had also not finished. | A finished install walk, with the details given to the walker directly. |
| C | The upkeep walk had no installation to look after, because the install walk had not finished. | A copy of a finished install walk's outputs. |
| W | A walk rule forbade sourcing settings files, and one quick-test step sources one. | Allow that one documented step in the walk rules. |
| B | The clean home folder put the system's old bash first, and the briefs said Java was present when it was not. | A clean home with the baseline a real user has. |

Three more limits apply to every row:

- The container images were already cached on the test machines. First-pull time and pull failures were not seen.
- The phone kit's planted control was not engaged. Both phone walkers went around it, so the phone kit has not yet proven it reports a stuck point. The desktop and shell controls were caught.
- Rounds 2 and 3 (bad input, failures, phone and agent parity) were planned but not walked.

## What held

- On a desktop or a phone, a new user learned the supported inputs, the software, the memory and disk needs and a small public test. It took under 20 actions, in both walks.
- From the docs, a new user found the demo and understood a finding. It was CYP2D6 \*1/\*3, intermediate metaboliser, with its evidence, its caveat and the advice to confirm it clinically. The demo says plainly that its people are invented.
- The demo showed ready, waiting, skipped and finished states during the replay.
- In the demo family, a sensitive result stayed closed until opened, closed again on request and stayed closed on return.
- A person traced a variant to its chromosome, band and position. A search for an absent gene said it covers only the reported results.
- On a desktop, the demo report downloaded and opened.
- The coding agent with the docs index planned a VCF route with its limits and outputs. The freshness command gave a clear message when `og` was missing.
- The verifier found the README badge and the public test data release working.

## What an agent cannot judge

- Feelings: trust, worry or relief when a person reads a health result about themselves.
- A real screen reader. The run used the accessibility tree and contrast checks as stand-ins.
- Real devices: touch on a real phone, saving or printing a report on a phone, other browsers.
- Numbers that need real users: how long a person takes to understand a result, and how many give up.
- Clinical meaning. Invented and partial data test the mechanics, not whether a negative result can be trusted.

## Walk it yourself

Pick a journey above, start from its entry point and try its goal on your own machine. Then open an issue at https://github.com/GeiserX/OwnGenome/issues.

Tell us where you paused, guessed, waited, or saw a result you did not trust.
