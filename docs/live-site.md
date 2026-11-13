# Watch a Run Live

A full run takes a day or more. The live site shows each person's results as their steps finish, in a browser on any device that can reach the machine. Each genome reads as a report: a lead paragraph that says what to take to a doctor, what to tell a prescriber and what is worth a look, then ten numbered sections from what matters most to what the run cannot see, each finding with a page of its own. While the run goes, a short band under the name counts the steps, and the timeline of the steps closes the page, folded. A home page lists everyone by family, and each family has a report of its own. The pages repaint as messages arrive; nobody has to reload them.

`bin/owngenome_live.py` is the server. It uses Python's standard library only, so the `python3` on the host is enough, and it has its own password, access log and limits, so it needs no web server in front of it.

## Try the demo

https://geiserx.github.io/OwnGenome/demo/ is the live site on DEMO-001, the invented sample of the docs pictures. It opens on the finished report of a made-up run and replays it from the first task to the database re-check after it: press Play, or drag the slider, and watch the cards fill in as their steps end. Every number, result and family member in it is made up. The family is an invented trio: DEMO-001, the partner DEMO-002 and their son DEMO-003. Each of the three has a finished run with its own findings, so every card widget and the whole family page have something to draw: DEMO-002 carries an FMR1 premutation-range count, a CHIP call and a MoChA loss; DEMO-003 a class 4 AnnotSV deletion and a Y path.

![The person page of DEMO-001 at the end of the demo run: the name with its family and run state, the contents of the ten sections with a count per level, the lead paragraph (no finding for a doctor, two medicine results, six worth a look), the five-level ledger, the data-quality line and the start of the chromosome map of the findings](images/demo-live-site.png){ loading=lazy }

The pages and their plugins are the files of `assets/live/`, unchanged. Each moment of the replay is the server's own `state.json` at that moment: `tests/demo/build_live_demo.py` loads every plugin, stores the made-up messages of `tests/fixtures/live/demo-001.weblog.jsonl` one by one as the ingest port would, and reads the state after each. The cards come from the `results.json` the server would give at the end of the run; each shows once the replay reaches the end of its step. The home page, the charts, the search box, the run-health card, the database re-check and the family page work as they do on your own site. The genome map of each person page is the document `og map build` writes for that run, built once when the demo is. The genome viewer does not work: it opens your own reads on your own computer, so its links lead to a page that says it is not in the demo.

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

With `OWNGENOME_LIVE_URL` set, `run-all.sh` writes `<sample>/nextflow/live.config` and passes it to Nextflow with `-c`. That config loads the nf-weblog plugin (pinned as `NF_WEBLOG_PLUGIN` in `versions.env`; Nextflow downloads it on the first live run) and points it at `OWNGENOME_LIVE_URL/ingest/<sample>`, so every task's progress is posted. The numbered script steps that run after the pipeline (GRIDSS, imputation prep, the somatic caller, the extra callers) post when they start and how they end (the benchmark, `benchmark-variants.sh`, has no step number and posts nothing), and `run-all.sh` posts once more as it ends, with its exit code. A site that is down or slow costs one warning; the run never waits for it and never fails for it. `OWNGENOME_LIVE_URL` must be a plain `http://host:port` URL. Set `OWNGENOME_LIVE_PAGE_URL` to the address you open the site at (`http://127.0.0.1:8787`, or the name behind your proxy) and `run-all.sh` also prints the run's page, `OWNGENOME_LIVE_PAGE_URL/runs/<sample>/`, as it starts; it prints it only together with `OWNGENOME_LIVE_URL`.

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

**Home** (`/`). A date line, then one sentence on where everyone stands (how many reports are finished or being written, and whether any finding needs a doctor), and beside it the findings that need attention across everyone, one line per person and level, each opening its finding page. A run still going comes first, as "Being written", with its steps done and the sections already written. A run with no end that has written nothing for a long time reads "Stopped" instead, with the day it was last heard from, and comes after the runs still going: twice the longest step's usual hours (27 h, as DeepVariant's 13.5 h can pass with no message) and never under 6 h. It reads "Being written" again as soon as something new arrives (a message, a trace row, a status row). Then each family: a small pedigree in the members' colours and one tile per person, with the first lines of their report's lead and their counts by level. A family member with no run on this server is listed, dashed. Last, every run as a log: state, steps, hours, when it ended. Without scripts the page is a plain list of runs.

**A person** (`/runs/<run>/`). The name, the sex the reads show against the declared one, the person's place in the family and the run state. Down the left, the contents: the ten sections with how many findings each level has, then the step timeline and the family. The page opens on the lead: a paragraph built from the findings (data-quality problems first, then the findings for a doctor named one by one, the medicine genes, the carrier genes, what is worth a look, and the mean depth; with no finding for a doctor it says so only when a clinical screen ran, and otherwise that the clinical screens have no result yet), the ledger of the five levels with a count each, the data-quality line, and Figure 1, the genome map. Then the sections:

1. **To act on**, 2. **Medicines**, 3. **For your family**, 4. **Worth a look**, 5. **Can I trust this?**: each finding of the level, with what it is, what to do, and in the margin its evidence, its caveat and its link to the report; then the cards of the report's sections that belong there.
6. **Origins** and 7. **Traits and scores**: the ancestry and trait cards.
8. **Assessed, no finding**: one line per section that ran and raised nothing, each opening its full card.
9. **What it cannot see**: what the run does not assess and the steps that did not run.
10. **Sources and files**: the database versions, the database re-check, and the report and its files.

A card has its rows (the first ones, then a button for all), its own sentences and its charts at reading width, and in a margin beside them its numbers, its caveat (the first sentence, then the rest behind "Full caveat") and its links; below 1100 px the margin goes under the text. A card with more than six lines keeps the ones that say something in view and folds the rest (an allele absent, a call not determined) behind one line that counts them. A section waiting for its step says "being written" or "waiting", one whose step failed says so, and a newer database release shows on the cards it would change. Some cards draw their numbers as well: the repeat expansions as one gauge per locus (the count against Stranger's normal and pathologic thresholds, with its 95% interval; an allele not sized from short reads is drawn hollow and never called normal), Lp(a) as the KIV-2 count with its 95% interval, CHIP as each call's VAF against the detection limits at this depth (hatched below every gene's limit, shaded where some genes see a clone and others do not, and each call named against its own gene's limit), the copy numbers of SMN1, SMN2, HBA1, HBA2 and the RCCX module as dots (Parascopy's SMN1 and SMN2 at exon 7 too, from short reads, hollow with its reason when the count is not resolved), and the Y haplogroup as its path from the root, with the marker counts of each level. Each drawing has the same numbers as text beside it. APOE stays behind its warning in the traits card; the polygenic scores show only once they are placed against the ancestry panel. A child's report has one gate at the data loader. While any result to act on is unopened, the page shows fixed explanatory text, the number of act results and one opaque opening control for each. An act opened on purpose shows only its own explanation, plan, evidence and caveat. It has no related list, next finding, map mark or chart. Nothing else from that child's report reaches a view: no cards, tables, level counts, images, file links, search box, search status or progress headlines.

The cost is that a child with an unopened act result sees none of the rest of the report until every act result is opened. Opening the last act releases the full report. Closing an act from the map replaces the document with the closed view. That replacement discards search results, view caches and in-flight view callbacks. A short-lived, one-use opening ticket contains only opening flags and a digest, and is consumed immediately. Reloading closes the acts again. Map and section addresses use the closed controls while the gate is shut. After the last opening, the page returns to its default view.

The controls use ordinal handles such as `closed-1` and `opened-1`. Guessed finding hashes use the same first closed handle. Opening or closing a control preserves the addressed act as the ordinals change. A finding address for an opened act shows that act alone while the gate is shut. Other finding addresses show the fixed explanation and controls. Adults and children with no act result use the full report. Unknown family roles and a child's unreadable results stay closed.

The home tile and attention row use fixed text for a child whose gate is shut. The family loader removes that child's report objects before rendering. Every retained side-by-side row has a fixed closed marker in that child's cell. Rows supplied only by that child are absent. Parents keep their values. Each child opens independently. The family tree remains visible because it comes from the declared pedigree, not the child's report.

This replaces the per-element veil. There is no fixed-panel list, gene-matching rule or prose exemption registry. `results.veil.acts` lists the act identities. The loader derives the controls from their own records. The map retains its separate APOE warning on every full report, including adults' reports.

This is a browser-only boundary. Static result, map, chart, search and family JSON files remain downloadable and contain their complete data. The gate protects the rendered document, attributes, accessibility text, title, address, canvas, SVG, browser storage and requests made by views. It is not an access-control boundary for those static files.

- **How the run went.** One line per step the run planned, on a time axis: a bar from its start to its end, or to now while it runs, with an outline to its usual hours. The colour is the step's state from `bin/run_state.py` ([Where a Run Stands](run-state.md)): waiting, running, done, failed, skipped, not finished. Under it the tasks running now against their usual time, a step's detail, the wire of what was found so far and the plugins' cards (run health, the database re-check).
- **A step's detail.** Click a step: its state and reason, start and end, the usual hours, each process with its tasks waiting, running, done, from the cache or failed, and the headline of its sections with a link into the report.
- **Found so far (the wire).** Newest first, one line per finished section that says something (a count, a warning, a failure), worded by the same code that writes the report and cut to its first sentence; the step's detail has the whole text. The sections whose step only finished fold into one line ("Variant calling and 21 more steps finished since 13:54"). The eight newest lines show, the rest sit under "Every event". The ClinVar line, with no hit, is information, never a clean result.
- **The report.** Once written, `report.html` opens the run's `<sample>_report.html`, and the files it links to beside it.
- **The genome viewer.** `/runs/<run>/viewer/` shows the run's reads and variants in igv.js, as [`scripts/owngenome-view.sh`](viewer.md) does on one computer. The report's Viewer links, the charts' "See the reads" links and a finding page's reads open it at their position.

The header search appears only on a full report. It finds a gene or rsID across the results and opens the card that names it. A result row with a field of the wrong type produces the finding "A row of this section could not be read". It takes the highest level its rule can give. If that level is act, it keeps a child's gate shut until someone opens it. On a full report, the card says which table cannot be drawn and the other cards remain available. `#s<n>` opens a section, `#s/<section>` a card (the report's own anchors), `#f/<finding>` a finding, `#q/<chapter>` the section a report chapter falls in, `#progress` the timeline.

**Figure 1, the genome map.** The 24 chromosomes as rows, to scale (23 for a person with no Y), each with its bands, the centromere as a waist and the stretches few reads can map hatched (centromeres, heterochromatin, the acrocentric arms). Each finding with a place sits on its chromosome as its level's shape, with its gene beside it when the name fits; runs of homozygosity of 5 Mb or more are a line under the row, shorter ones a faint band. A finding that names a gene but no position (a medicine gene, Lp(a)) sits on that gene's span; a finding with neither stays in its section only. Click a chromosome (or, with the map focused, press `[` or `]`) and it opens in place: an overview with the part in view, a ruler in Mb, the bands with their names, a findings lane with each gene and its level's words, and a lane for the runs of homozygosity, each run with its length. Drag to move along it, scroll or press the + and − buttons beside the zoom stops to zoom in (as far as about 10 Mb across a wide screen; closer views come later), and press Esc (or click Genome) to go back. Back in the browser returns to the view before, the whole genome when the address had no place. A dashed line marks the centre, and the margin says what is under it: the chromosome, the position and the band, the findings in view and the runs. The findings count, the chromosome sentences, the runs sentence and each finding's phrase come from the map document; the page adds its own short lines: a lane or stretch with no finding, the runs in view, how to move, and why a layer is not on the map. A replay draws the runs of homozygosity, and says the sentences about them, only once their step has ended. The address follows the map: `#genome`, `#chr7`, or `#chr7:117.4-137.4M`, and a link to one of them opens that view; a place the page cannot read shows a short note and leaves the map as it was. The map comes from `/runs/<run>/map.json`, which the map plugin (`bin/owngenome_live_map.py`) builds with `og map build` from the run folder and the results the page gets (`rust/crates/og-map/README.md` describes the document). It needs og on the host (`scripts/setup.sh --og`); without it, or in a browser that cannot unpack the file (no `DecompressionStream`), the figure says so and lists the findings with a place in a table instead. A map that fails to load again after it was drawn says so the same way and draws nothing more from the map it had. A layer whose file is missing or cannot be read says so in words, never as a blank lane.

**A finding** (`/runs/<run>/finding/<id>/`). One finding on a page of its own: its level and how many others share it, the title, the gene and variant, a plain explanation, and **What to do** as numbered steps. Then the evidence (its numbers, how common the change is on a log scale from one in a million to everyone, the caveat in the margin), the reads with a link to the genome viewer, the rest of the family, the related findings of the same gene or step, and the section of the report it comes from. The contents on the left lead back to the report and on to the next finding. The page is 404 unless the id is a finding of that run's `results.json` now.

**A family** (`/runs/<family>/`). Five chapters show the pedigree, measured relationships, couple carrier screens, each child's inheritance and the members side by side. The shared data gate above applies before the family view receives child report data. Gated children have fixed opening controls and closed cell markers. Parent values and independently opened children remain available. The full family report returns for a child after that child's final act opens.

On a phone each page is one column and each table row becomes a short list of label and value. On a person's report the contents give way to a tab bar at the bottom (To act, Medicines, Family, Trust, More); under 700 px the report is short: the findings of sections 1 to 4 (to act on, medicines, family, worth a look) stay open, while the cards behind them, every information card and the bodies of sections 5 to 10 start closed, each behind one line that says what it holds (section 5 stays open when it holds a finding about data quality, or a card the report flags with no finding of its own). The report gate above applies at every width. A fold you open stays open when the page is drawn again, the phone's folds too when you turn the phone to a wide screen and back. The counts read one line per level, and the "How to read this" note and the links between findings are left out. A link to a card, a section or a finding, a search hit and a tab open the folds around it (a Ctrl, Cmd or Shift click leaves this page as it is); a finding page drops its contents (the header links back to the report), and a family page keeps its chapter list above the chapters. The pages are dark. `tests/demo/render_pictures.py --out DIR --only live-phone` draws the picture below from the demo.

![The person page of DEMO-002 on a 390 px phone: the lead, the counts one line per level, the two medicine results and the nine results worth a look open, every other card one line with a plus, sections 5 to 10 each folded to one line that names what they hold, and the tab bar at the bottom](images/demo-live-phone.png){ loading=lazy width="390" }

## The genome viewer

`bin/owngenome_live_viewer.py` serves it, behind the same password as every page. `owngenome-live.sh` passes it the reference (`REF_FASTA`, as in every step) and igv.js, which the first start downloads into `GENOME_DIR/tools/` and checks against `IGV_JS_SHA256` in `versions.env`. Without the reference's `.fai` or igv.js the site starts without the viewer and says why; `OWNGENOME_LIVE_VIEWER=off` leaves it out on purpose. Running `owngenome_live.py` yourself (in Docker, say), pass `--viewer-fasta` and `--viewer-igv-js`.

Per run it serves the folder's BAM or CRAM, its VCF and its SV VCF, each with its index, and the reference: the files `owngenome-view.sh` picks by default. A track whose file or index is a symlink out of `GENOME_DIR` is left out. The page and the files carry the viewer's Content-Security-Policy, so the page talks to this site only, and igv.js's calls to outside services are blocked.

The reads cross the network to the browser that logged in: only the byte ranges igv.js asks for, but they are the sample's reads. Off loopback, use TLS (`--tls-cert` and `--tls-key`, or a proxy that ends TLS); without it the password and the reads cross the network in clear. Opening the viewer writes one access-log row; the byte ranges igv.js reads write none.

`GET /` lists the runs. The home page reads `/runs.json`, then resolves family membership and each child's acts before projecting the run summaries. A person's loader reads `results.json` and the family routing information. Full reports contain the chapters, findings, journal, lead and source metadata built by `bin/live_results.py`. The gate passes none of that report to views while an act remains closed, apart from explicitly opened acts' own records. The family loader reads `family.json` and each child's results, and refreshes every 30 seconds. A changed act record closes again. An unchanged record keeps its opening for the current page.

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

The page has a replay mode for a demo with no server: `<html data-mode="replay">` reads `frames.json` (`{"interval_ms": 1000, "frames": [{"t": ..., "state": <state.json>}]}`) and plays it with a play button and a slider, through the same drawing code. A replay whose last frame is a finished run opens at that frame, on the finished report, until you press Play, which starts again from the first frame; with `"autoplay": true` it plays from the start at once. The bar shows each frame's moment as day, month and time in UTC ("2 Jan 10:05"), the full time on hover.

`tests/test_owngenome_live.py` checks the server, login, message storage, plugins and browser. `tests/test_live_results.py` checks card caveats, wording and login boundaries. `tests/test_live_journal.py` checks finding levels, leads and routes. `tests/test_live_demo.py` checks the demo, maps, phone layouts and report gate. `tests/test_owngenome_live_family.py` runs family checks and its browser property in a separate CI shard. `tests/test_live_post.sh` and `tests/fake-docker/run-all-live.sh` check progress messages from the launcher.
