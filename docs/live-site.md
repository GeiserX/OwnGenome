# Watch a Run Live

A full run takes a day or more. The live site shows each person's results as their steps finish, in a browser on any device that can reach the machine. Each genome reads as a report: a lead paragraph that says what to take to a doctor, what to tell a prescriber and what is worth a look, then ten numbered sections from what matters most to what the run cannot see, each finding with a page of its own. While the run goes, a short band under the name counts the steps, and the timeline of the steps closes the page, folded. A home page lists everyone by family, and each family has a report of its own. The pages repaint as messages arrive; nobody has to reload them.

`bin/owngenome_live.py` is the server. It uses Python's standard library only, so the `python3` on the host is enough, and it has its own password, access log and limits, so it needs no web server in front of it.

## Try the demo

https://geiserx.github.io/OwnGenome/demo/ is the live site on DEMO-001, the invented sample of the docs pictures. It replays a made-up run from the first task to the database re-check after it: press Play, or drag the slider, and watch the cards fill in as their steps end. Every number, result and family member in it is made up. The family is an invented trio: DEMO-001, the partner DEMO-002 and their son DEMO-003. Each of the three has a finished run with its own findings, so every card widget and the whole family page have something to draw: DEMO-002 carries an FMR1 premutation-range count, a CHIP call and a MoChA loss; DEMO-003 a class 4 AnnotSV deletion and a Y path.

![The person page of DEMO-001 at the end of the demo run: the name with its family and run state, the contents of the ten sections with a count per level, the lead paragraph (no finding for a doctor, two medicine results, six worth a look), the five-level ledger, the data-quality line and the start of the chromosome map of the findings](images/demo-live-site.png){ loading=lazy }

The pages and their plugins are the files of `assets/live/`, unchanged. Each moment of the replay is the server's own `state.json` at that moment: `tests/demo/build_live_demo.py` loads every plugin, stores the made-up messages of `tests/fixtures/live/demo-001.weblog.jsonl` one by one as the ingest port would, and reads the state after each. The cards come from the `results.json` the server would give at the end of the run; each shows once the replay reaches the end of its step. The home page, the charts, the search box, the run-health card, the database re-check and the family page work as they do on your own site. The genome viewer does not: it opens your own reads on your own computer, so its links lead to a page that says it is not in the demo.

What differs from your own site:

- It is a static folder published with these docs, so it has no password and no access log. Nothing in it is private.
- The run-health card measures an invented machine (64 GiB of RAM, a 2 TiB disk 55% used), not the one that built the page.

## Start it

```bash
export GENOME_DIR=/path/to/your/data
./scripts/owngenome-live.sh hash-password                         # once: asks a password, stores only its hash
./scripts/owngenome-live.sh                                       # the site, on http://127.0.0.1:8787
OWNGENOME_LIVE_URL=http://127.0.0.1:8788 ./scripts/run-all.sh your_name male   # in another terminal
```

`./scripts/owngenome-live.sh your_name` does the same and prints the address of that run's page. `hash-password` writes `GENOME_DIR/.owngenome-live/password` with mode 600. The file holds a PBKDF2-SHA256 hash, never the password. Open http://127.0.0.1:8787, log in, and pick the run.

With `OWNGENOME_LIVE_URL` set, `run-all.sh` writes `<sample>/nextflow/live.config` and passes it to Nextflow with `-c`. That config loads the nf-weblog plugin (pinned as `NF_WEBLOG_PLUGIN` in `versions.env`; Nextflow downloads it on the first live run) and points it at `OWNGENOME_LIVE_URL/ingest/<sample>`, so every task's progress is posted. The numbered script steps that run after the pipeline (GRIDSS, imputation prep, the somatic caller, the extra callers) post when they start and how they end (the benchmark, `benchmark-variants.sh`, has no step number and posts nothing), and `run-all.sh` posts once more as it ends, with its exit code. A site that is down or slow costs one warning; the run never waits for it and never fails for it. `OWNGENOME_LIVE_URL` must be a plain `http://host:port` URL.

A bare `nextflow run main.nf` can report too. Write a config like this one, with the sample id at the end of the URL (the page counts only the tasks tagged with it), and pass it with `-c`:

```groovy
plugins {
    id 'nf-weblog@1.2.0'
}
weblog {
    enabled = true
    url     = 'http://127.0.0.1:8788/ingest/your_name'
}
```

## What the pages show

Every result carries one of six levels, with a word, a shape and a colour, so the shape alone tells them apart: **Take to a doctor** (◆), **Tell your prescriber** (a pill), **Worth a look** (▲), **For family planning** (two rings), **Data quality** (✓, red when the data has a problem) and **For information** (○). `bin/live_results.py` sets the level from the rows the report already prints: an ACMG secondary finding, a pathogenic CPSR row or two hits in a carrier-list gene is to act on; a ClinVar hit, an AnnotSV class 4, a repeat outside Stranger's normal, a CHIP or MoChA call or an Exomiser top gene is worth a look (a full-mutation repeat, an AnnotSV class 5 and a positive VNtyper call are to act on); a carrier-list row, a gene both partners carry and one SMN1 copy are for family planning; a CPIC gene with changed function is for the prescriber; a sex mismatch, a contamination warning or a second run of the same person is data quality. Everything else is for information. No page says a section is fine: a section with no finding says "No finding in what these steps assessed" and what it looked at.

**Home** (`/`). A date line, then one sentence on where everyone stands (how many reports are finished or being written, and whether any finding needs a doctor), and beside it the findings that need attention across everyone, one line per person and level, each opening its finding page. A run still going comes first, as "Being written", with its steps done and the sections already written. Then each family: a small pedigree in the members' colours and one tile per person, with the first lines of their report's lead and their counts by level. A family member with no run on this server is listed, dashed. Last, every run as a log: state, steps, hours, when it ended. Without scripts the page is a plain list of runs.

**A person** (`/runs/<run>/`). The name, the sex the reads show against the declared one, the person's place in the family and the run state. Down the left, the contents: the ten sections with how many findings each level has, then the step timeline and the family. The page opens on the lead: a paragraph built from the findings (data-quality problems first, then the findings for a doctor named one by one, the medicine genes, the carrier genes, what is worth a look, and the mean depth; with no finding for a doctor it says so only when a clinical screen ran, and otherwise that the clinical screens have no result yet), the ledger of the five levels with a count each, the data-quality line, and Figure 1, the chromosomes with each finding at its place. Then the sections:

1. **To act on**, 2. **Medicines**, 3. **For your family**, 4. **Worth a look**, 5. **Can I trust this?**: each finding of the level, with what it is, what to do, and in the margin its evidence, its caveat and its link to the report; then the cards of the report's sections that belong there.
6. **Origins** and 7. **Traits and scores**: the ancestry and trait cards.
8. **Assessed, no finding**: one line per section that ran and raised nothing, each opening its full card.
9. **What it cannot see**: what the run does not assess and the steps that did not run.
10. **Sources and files**: the database versions, the database re-check, and the report and its files.

A card has its rows (the first ones, then a button for all), its own sentences and its charts at reading width, and in a margin beside them its numbers, its caveat (the first sentence, then the rest behind "Full caveat") and its links; below 1100 px the margin goes under the text. A card with more than six lines keeps the ones that say something in view and folds the rest (an allele absent, a call not determined) behind one line that counts them. A section waiting for its step says "being written" or "waiting", one whose step failed says so, and a newer database release shows on the cards it would change. Some cards draw their numbers as well: the repeat expansions as one gauge per locus (the count against Stranger's normal and pathologic thresholds, with its 95% interval; an allele not sized from short reads is drawn hollow and never called normal), Lp(a) as the KIV-2 count with its 95% interval, CHIP as each call's VAF against its gene's detection limit, the copy numbers of SMN1, SMN2, HBA1, HBA2 and the RCCX module as dots, and the Y haplogroup as its path from the root, with the marker counts of each level. Each drawing has the same numbers as text beside it. APOE stays behind its warning in the traits card; the polygenic scores show only once they are placed against the ancestry panel. On a child's page in a family, a result to act on stays closed behind "Show the result" until someone opens it on purpose: its title, gene and variant stay out of the lead, the ledger, the chromosome map, the home page and its own page, and the card it comes from stays closed with it (its rows, numbers and drawing), in its section and in section 8. Each result opens on its own: opening one leaves every other result to act on closed, on the report and on the finding's page. A person in a family whose family report has not answered yet, or answered without listing them, is treated the same way, since the page cannot tell whether they are a child.

While the run goes, the page still opens on the report being written. A short band under the name says how far the run is ("23 of 27 are done, about 8–13 hours remain"), with the steps done as a number and a bar, the time elapsed, the time to go from the usual step times and the retries, and links to the step timeline. Each section and each line of the contents says whether it is written, being written or waiting. The run itself, below, closes the page after section 10, folded until opened (or reached through `#progress`, the band or the contents); a page with no results yet keeps it open.

- **How the run went.** One line per step the run planned, on a time axis: a bar from its start to its end, or to now while it runs, with an outline to its usual hours. The colour is the step's state from `bin/run_state.py` ([Where a Run Stands](run-state.md)): waiting, running, done, failed, skipped, not finished. Under it the tasks running now against their usual time, a step's detail, the wire of what was found so far and the plugins' cards (run health, the database re-check).
- **A step's detail.** Click a step: its state and reason, start and end, the usual hours, each process with its tasks waiting, running, done, from the cache or failed, and the headline of its sections with a link into the report.
- **Found so far (the wire).** Newest first, one line per finished section that says something (a count, a warning, a failure), worded by the same code that writes the report and cut to its first sentence; the step's detail has the whole text. The sections whose step only finished fold into one line ("Variant calling and 21 more steps finished since 13:54"). The eight newest lines show, the rest sit under "Every event". The ClinVar line, with no hit, is information, never a clean result.
- **The report.** Once written, `report.html` opens the run's `<sample>_report.html`, and the files it links to beside it.
- **The genome viewer.** `/runs/<run>/viewer/` shows the run's reads and variants in igv.js, as [`scripts/owngenome-view.sh`](viewer.md) does on one computer. The report's Viewer links, the charts' "See the reads" links and a finding page's reads open it at their position.

The header search finds a gene or rsID across the results and opens the card that names it. `#s<n>` opens a section, `#s/<section>` a card (the report's own anchors), `#f/<finding>` a finding, `#q/<chapter>` the section a report chapter falls in, `#progress` the timeline.

**A finding** (`/runs/<run>/finding/<id>/`). One finding on a page of its own: its level and how many others share it, the title, the gene and variant, a plain explanation, and **What to do** as numbered steps. Then the evidence (its numbers, how common the change is on a log scale from one in a million to everyone, the caveat in the margin), the reads with a link to the genome viewer, the rest of the family, the related findings of the same gene or step, and the section of the report it comes from. The contents on the left lead back to the report and on to the next finding. The page is 404 unless the id is a finding of that run's `results.json` now.

**A family** (`/runs/<family>/`). A report in five chapters, opened by a lead: how many genomes are read together, whether the couple shares a carrier gene, whether the measured relationships match the tree, and how many children. I, the pedigree, drawn from `members.tsv` or `pedigree.tsv` in the members' colours, where each sequenced person with a run here opens their page and shows how many cards their report flags, beside the members and the analyses step 38 planned. II, who is related how: step 38's relationship check beside the shared segments, as a grid of every pair, with ⚑ where the measured relationship differs from the declared one. III, each couple's carrier screen in full (genes where both carry a finding, drawn as four children with the chance for each, X-linked genes, genes where one partner does), with the screen's caveat and a link to the couple report. IV, each child: the stretches inherited from each parent painted per chromosome (coloured by grandparent where step 38 could tell), the de novo candidates and the compound-het pairs, in trans first. V, the members side by side: their CYP2D6 calls, CPIC phenotypes, HLA types, haplogroups and blood groups, one column each. On a phone the grid becomes a list of pairs.

On a phone each page is one column and each table row becomes a short list of label and value. On a person's report the contents give way to a tab bar at the bottom (To act, Medicines, Family, Trust, More); under 700 px the report is short: the findings of sections 1 to 4 (to act on, medicines, family, worth a look) stay open, while the cards behind them, every information card and the bodies of sections 5 to 10 start closed, each behind one line that says what it holds. The counts read one line per level, and the "How to read this" note and the links between findings are left out. A link to a card, a section or a finding, a search hit and a tab open the folds around it; a finding page drops its contents (the header links back to the report), and a family page keeps its chapter list above the chapters. The pages are dark.

![The person page of DEMO-002 on a 390 px phone: the lead, the counts one line per level, the two medicine results and the nine results worth a look open, every other card one line with a plus, sections 5 to 10 each folded to one line that names what they hold, and the tab bar at the bottom](images/demo-live-phone.png){ loading=lazy width="390" }

## The genome viewer

`bin/owngenome_live_viewer.py` serves it, behind the same password as every page. `owngenome-live.sh` passes it the reference (`REF_FASTA`, as in every step) and igv.js, which the first start downloads into `GENOME_DIR/tools/` and checks against `IGV_JS_SHA256` in `versions.env`. Without the reference's `.fai` or igv.js the site starts without the viewer and says why; `OWNGENOME_LIVE_VIEWER=off` leaves it out on purpose. Running `owngenome_live.py` yourself (in Docker, say), pass `--viewer-fasta` and `--viewer-igv-js`.

Per run it serves the folder's BAM or CRAM, its VCF and its SV VCF, each with its index, and the reference: the files `owngenome-view.sh` picks by default. A track whose file or index is a symlink out of `GENOME_DIR` is left out. The page and the files carry the viewer's Content-Security-Policy, so the page talks to this site only, and igv.js's calls to outside services are blocked.

The reads cross the network to the browser that logged in: only the byte ranges igv.js asks for, but they are the sample's reads. Off loopback, use TLS (`--tls-cert` and `--tls-key`, or a proxy that ends TLS); without it the password and the reads cross the network in clear. Opening the viewer writes one access-log row; the byte ranges igv.js reads write none.

`GET /` lists the runs: each folder under `GENOME_DIR` holding `logs/run_status.tsv`, `nextflow/` or `run_manifest.tsv`, and each run the site has messages for. The home page reads `/runs.json`: each run's progress and attention count, and each family with its members. A run with results also has `levels` (its findings counted by level), `top` (up to four findings to show first, data-quality problems first: id, level, title), `lead` (its lead paragraph as plain text) `lead_closed` (the same with the names of the findings to act on left out, for a child's tile) and `screened` (whether ClinVar, the clinical filter or CPSR has a result; without one the home page never says that no finding needs a doctor). While a run goes these come from the summary its page last read, as long as no step has finished and `summary.json` has not changed since that read, so the home page never shows less than the person page. While a family's members cannot be read, the home page keeps the findings to act on closed for every run that no readable family lists.

The person page reads `/runs/<run>/results.json` from `bin/owngenome_live_results.py`: `bin/live_results.py` turns the run's summary (its `summary.json` when step 24 wrote it after every finished step, else the step files read again when a step ends; a rewritten `summary.json` is read again too, and the page fetches the cards again) into the report's chapters, one card per section. Each card also has its `level`; the file adds `findings` (each with its id, level, section, title, plain text, gene, evidence, caveat, plan and related findings), `levels` (the counts), `journal` (the ten sections and the cards in each), `lead` (the lead paragraph as text pieces, the bold ones marked, and the names of the findings to act on marked `act`), `level_words` and `databases`. It computes nothing about genetics and holds only what the report prints. The family page reads `/runs/<family>/family.json` from the family plugin, and asks again every 30 seconds, so a member's progress and a new couple screen show without a reload.

## Privacy and what is stored

- The site serves only its own pages, the run's state, the results the report prints (`results.json`, never `summary.json` itself), the report, the files the report links to by a relative link, and the genome viewer's files of the run (see above). Any other file is 404, and so is a link that leaves the run folder or a symlink out of `GENOME_DIR`. It never lists a folder: a link to a folder that holds one HTML file (the report's `cpsr/`, with CPSR's own report) opens that file, and any other folder link is 404.
- Each message is cut to what a progress view needs before it is stored: the task, its process and status, times, CPU and memory figures. Nextflow's parameters, command line and task scripts are dropped. Three paths stay on the server for later views (the output folder, the CPU and memory caps, the work directory); none of them reaches a page.
- Messages go to `<state dir>/runs/<run>/events.jsonl` (state dir: `GENOME_DIR/.owngenome-live`, or `OWNGENOME_LIVE_STATE_DIR`). The site writes nowhere else, so `GENOME_DIR` can be mounted read-only.
- The access log, `<state dir>/access.log`, has one tab-separated row per page or file read and per login (the genome viewer's files excepted: one row per page open, see above): time, client, the first 8 characters of a hash of the session, action, run, file, status. It never holds a query string, a request body, a result or the password. It moves to `access.log.1` past `--access-log-max-mb` (10).
- Every page carries a Content-Security-Policy that lets it talk to this server only, and the page loads no script from anywhere else. The run's own HTML (the report, CPSR, PharmCAT, MultiQC) is sent sandboxed: no script runs and nothing outside is fetched. CPSR's report draws its variant tables with scripts, so on the live site those tables stay empty: open the file from `<sample>/cpsr/` on disk to see them.

## The password and the limits

| Flag | Default | What it does |
|---|---|---|
| `--password-file FILE` | `OWNGENOME_LIVE_PASSWORD_HASH` | The hash to check logins against. Without one, the site listens on loopback only and asks for nothing; any other `--bind` exits 2 |
| `--bind`, `--port` | `127.0.0.1`, `8787` | The browser's address |
| `--ingest-bind`, `--ingest-port` | `127.0.0.1`, `8788` | Where Nextflow and `run-all.sh` post |
| `--session-hours` | 12 | A login lasts at most this long |
| `--session-idle-minutes` | 120 | A login unused this long ends; an open page counts as use only when it reads the state after a message, and its stream closes when the login ends |
| `--login-max-failures` | 5 | Failed logins per client in 15 minutes; then 429 until the window passes |
| `--max-streams` | 32 | Open live pages at once |
| `--max-body-kb` | 1024 | Largest message taken (413 above it) |
| `--allowed-host NAME` | none | A host name the browser may use besides loopback and the bind address (repeatable); any other Host header is 403 |
| `--trust-proxy CIDR` | none | Read `X-Forwarded-For` only from this proxy (repeatable) |
| `--tls-cert`, `--tls-key` | none | Serve HTTPS; the cookie is then `Secure` |
| `--cookie-secure` | off | Mark the cookie `Secure` when TLS ends at a proxy |
| `--viewer-fasta`, `--viewer-igv-js` | none | The genome viewer's reference and igv.js; `owngenome-live.sh` passes both |
| `--config FILE` | none | `key=value` lines setting any flag above (`session_hours=8`), read before the flags |

The session cookie is random (32 bytes), `HttpOnly` and `SameSite=Strict`. A POST from another site is refused.

The ingest port has no password: it takes messages from the processes of the run's own machine. That is why it listens on loopback by default: a local user can post fake progress there, but the port serves nothing, so nobody can read a run through it. A browser's post to it (one with an `Origin` header) is refused.

## Deploy it on a server

The same command runs on any Linux, macOS or WSL2 machine. Pick the reach you need.

**Loopback and ssh.** Leave the defaults and forward the port from your own machine:

```bash
ssh -L 8787:127.0.0.1:8787 user@your-server     # then open http://127.0.0.1:8787
```

**On the local network.** Set the password, then bind the machine's address:

```bash
OWNGENOME_LIVE_BIND=0.0.0.0 OWNGENOME_LIVE_ALLOWED_HOSTS="genome.example.org" ./scripts/owngenome-live.sh
```

`OWNGENOME_LIVE_ALLOWED_HOSTS` lists every name or address you will type in the browser: with `0.0.0.0` the site knows none of the machine's addresses, so `http://192.0.2.10:8787` is 403 until `192.0.2.10` is listed. Without TLS the password crosses the network in clear; use this on a network you trust, or add `--tls-cert` and `--tls-key`.

**In Docker.** The image is `PYTHON_IMAGE` from `versions.env`, the data read-only, the state folder writable outside it:

```bash
. ./versions.env
docker run --rm --user "$(id -u):$(id -g)" \
  -p 127.0.0.1:8787:8787 -p 127.0.0.1:8788:8788 \
  -v "$PWD":/pgp:ro -v "$GENOME_DIR":/genome:ro -v "$HOME/.owngenome-live":/state \
  "$PYTHON_IMAGE" python3 /pgp/bin/owngenome_live.py --root /genome --state-dir /state \
  --bind 0.0.0.0 --ingest-bind 0.0.0.0 --password-file /state/password
```

Inside the container both listeners bind every address, and `-p 127.0.0.1:...` keeps them on the host's loopback. Write the hash first: `mkdir -p -m 700 ~/.owngenome-live`, then `(umask 077; python3 bin/owngenome_live.py hash-password > ~/.owngenome-live/password)`, so the file is mode 600 from the moment it exists.

**Behind a reverse proxy.** Point the proxy at `127.0.0.1:8787`, pass the Host header, and turn off response buffering for `/runs/<run>/events` (the server sends `X-Accel-Buffering: no`, which nginx honours). Then:

```bash
./scripts/owngenome-live.sh --allowed-host live.example.org --trust-proxy 127.0.0.1 --cookie-secure
```

**As a service.** A systemd unit for the loopback setup:

```ini
[Unit]
Description=OwnGenome live site
After=network.target

[Service]
User=genome
Environment=GENOME_DIR=/path/to/your/data
ExecStart=/path/to/OwnGenome/scripts/owngenome-live.sh
Restart=on-failure

[Install]
WantedBy=multi-user.target
```

## Plugins

The core files (`bin/owngenome_live.py`, `assets/live/index.html`, `home.html`, `family.html`, `core.js`, `live.js`, `person.js`, `home.js`, `family-page.js`, `live.css`, `login.html`) are meant to stay small. A new view adds two files instead of editing them:

- `bin/owngenome_live_<name>.py` with `register(app)`. At start the server imports every such file, sorted, and calls it. `app.route(method, pattern, fn, auth=True)` adds a page or an endpoint (behind the login unless `auth=False`), `app.run_file(glob)` lets the browser fetch more files of a run folder, `app.report_filter(fn)` rewrites the served report, `app.state_hook(name, fn)` adds `plugins[name]` to `state.json`, `app.run_kind(fn)` lists more runs, `app.add_argument(...)` adds a flag, `app.log(...)` writes an access-log row, and `app.config` holds the parsed flags.
- `assets/live/plugins/<name>.js`. The run page loads every script in that folder after `live.js`; it calls `window.PGP.register({name, onState, onEvent, panel})`. `onState(state)` gets each new `state.json`, `onEvent(message)` each message from the stream, and `panel(element)` gets a card of its own under the run timeline, or in the element whose id `slot` names (the search box sits in the header this way). A family's page loads the same scripts without `live.js`, so a plugin returns early when `window.PGP` is missing.

The page has a replay mode for a demo with no server: `<html data-mode="replay">` reads `frames.json` (`{"interval_ms": 1000, "frames": [{"t": ..., "state": <state.json>}]}`) and plays it with a play button and a slider, through the same drawing code.

`tests/test_owngenome_live.py` checks the server, the login, the stored messages, the plugin interface and the page in a headless Chrome; `tests/test_live_results.py` checks `results.json` (every card's caveat word for word, no wording that reads as an all-clear, `results.json` behind the login); `tests/test_live_journal.py` checks the levels, the lead paragraph and the finding page's route; `tests/test_live_demo.py` builds the demo and checks every page, phone width included; `tests/test_live_post.sh` and `tests/fake-docker/run-all-live.sh` check what `run-all.sh` posts.
