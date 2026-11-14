# The Genome Map

The genome map is Figure 1 of each person's page on the live site (`docs/live-site.md`). It draws one person's genome from the whole set of chromosomes down to the exons of one gene, with the findings, the runs of homozygosity and the genes in place. It runs in the browser from one file the server builds for the run, and it needs no outside service: gene names and gene models come from the GENCODE files setup installs.

## What the map shows

The map has five levels. The level follows how many bases one screen pixel covers, so each level shows the same detail at any window size; the reads level is any view under 50 kb.

| Level | When | What it draws |
|---|---|---|
| Genome | the whole genome | the 24 chromosomes as rows, to scale (23 for a person with no Y), with their bands, the centromere as a waist and the stretches few reads can map hatched; each finding with a place as its level's shape, with its gene beside it when the name fits; runs of homozygosity of 5 Mb or more as a line, shorter ones as a faint band |
| Chromosome | more than 16,667 bases per pixel | one chromosome open in place: an overview with the part in view, a ruler, the bands with their names, a findings lane with each gene and its level's words, and the runs of homozygosity with their lengths |
| Region | 833 to 16,667 bases per pixel | the same lanes, with each band named with its chromosome (`16p13.3`), and the genes lane: every gene that fits in three rows, its name beside it; the genes a finding names come first, in bold, then the protein-coding genes, then the rest |
| Gene | 833 bases per pixel or fewer, a view of 50 kb or more | the genes as models: the span as a thin line with its strand's chevrons and the exon blocks inside it; each finding is a dashed tick down through every lane, with its change (HGVS) under its gene |
| Reads | a view under 50 kb | drawn as the gene level |

A dashed line marks the centre of an open chromosome. The margin beside the map says what is under it: the chromosome, the position and the band, the findings in view and the runs of homozygosity.

The exon blocks are the union of the exons of every basic GENCODE transcript of the region: a block means "an exon of some transcript here", not the exons of one transcript.

## Moving around

- Click a chromosome to open it, and drag to move along it. Scroll, the + and − buttons or the + and − keys zoom.
- The zoom stops above the map (Genome, Chromosome, Region, Gene) jump to a level and keep the place in view. From the whole genome they go to the first finding by level.
- The path under the title (`Genome › chr16 › 16p13.3 › 2.050-2.300M`) shows where the map is; each step up is a button.
- Click a gene to zoom to it. Click a finding to open its card in the margin: its level, title, gene, position, change and the map's sentence about it, with a link to the finding's own page.
- When an open chromosome has no finding in view, a "Nearest" button names the closest finding on that chromosome, or the first one on the next chromosome that has one, and flies there.
- The address follows the map: `#genome`, `#chr7`, `#chr7:117.4-137.4M`, `#chr7:117.480-117.669M` or `#chr7:117559001-117599000` (1-based, both ends included). A link to one of them opens that view, and Back returns to the view before.

### Keys

| Key | What it does |
|---|---|
| `/` | puts the cursor in the search box |
| ← → | move a tenth of the view; with Shift, a whole view |
| + − | zoom in or out 2×, at the centre |
| `[` `]` | the chromosome before or after |
| `0` | the whole genome |
| Esc | closes the card, then steps out one level: reads to gene, gene to region, region to the whole chromosome, chromosome to the genome |
| Tab | steps through the findings in view, each with a ring on the map and its card, then leaves the map |
| `L` | opens the layers and moves to them |
| `T` | folds or opens the In view table |

The keys work with the map focused (`/`, `L` and `T` anywhere on the page) and never inside a text box or with Cmd, Ctrl or Alt.

### The whole-genome strip

While a chromosome is open, a strip above the map draws every chromosome end to end, 5 Mb apart, with a tick for each finding on the map. The view is a brass brush, at least 3 pixels wide so a view of a few bases still shows, and the rest of the strip is dimmed. Click a chromosome on the strip to open it whole. Drag across the strip to open the stretch you brush, on the chromosome where the drag started; its size shows beside the brush while you drag. Drag the brush itself to move along the open chromosome. On a phone the strip is a picture only.

### The In view table

Under the map, the "In view" table lists every mark the map draws one by one in the view: findings, runs of homozygosity and, from the region level, genes. Its columns are Kind, What, Where (1-based), Size, Who, Copy and Source. It starts sorted by place with the findings first; click a column header to sort by it, and again to reverse. It shows 12 rows and a "Show all" button, which lists every row (a scrolled list past 200 rows, at most 5,000; past that the table says to zoom in). A layer the map draws as a density at this level gets one line under the table instead of rows, such as the shorter runs of homozygosity at the genome level. Click a row, or its name, to fly the map to it; a finding's card opens too.

The table follows the map once it stops moving, and tells a screen reader only how many marks are in view, once each time that number changes. On a phone it starts folded under its count. The canvas, the strip and the table use one rule for what the map draws one by one at each level (`drawsAt()` in `assets/live/map.js`), and the live demo's tests hold it to the same rule in `og map build --predicate-table`.

### Search

The search box finds a gene by its name (every GENCODE basic gene, not only the protein-coding ones), a band by its name (`13q13`), a finding drawn on the map by its gene, its id or words of its title, and a place written as the address writes it. The gene names are fetched with the first key typed. Pick a line with the arrow keys and Enter, or click it: the map flies there. Search offers only what the map draws, so a result the page keeps closed is never suggested.

## Layers

| Layer | From | Notes |
|---|---|---|
| Bands | `assets/reference/cytoBand.hg38.txt` | the UCSC GRCh38 cytoband table the repository ships |
| Findings | the results the page shows | placed from each finding's position; a finding that names only a gene sits on that gene's span; a finding with neither stays in its section |
| Runs of homozygosity | step 11 (`<run>/vcf/<run>_roh.txt`) | runs of 5 Mb or more drawn solid; a run mostly inside a hatched stretch, or on a male's X or Y, is marked as not a real run |
| Genes | the GENCODE files of setup | names and models from the region level down; without the files the lane says "gene names need step 08's GENCODE file" |

A layer whose file is missing or cannot be read says so in words in its lane and in the layer list, never as a blank lane. "Paint the copies" stays at None on a person's page: which parent each stretch came from needs a family run, and the control says why it is off.

## Closed results

On a child's page in a family, a result to act on stays closed until someone opens it on purpose, and APOE stays behind its warning on every page. While a result is closed, the map draws no mark for it, the strip has no tick for it, the In view table has no row for it, its card cannot open, search does not suggest it and Nearest does not lead to it. A dashed chip above the map counts the closed results without naming a gene or a place. `docs/live-site.md` describes the rule.

## The files behind the map

The server builds the map's files with `og map`, one GET request each, behind the site's login:

| Route | What it is |
|---|---|
| `/runs/<run>/map.json` | the person's map document: `og map build` with the results the page gets |
| `/map/genes/<chrom>.json` | one chromosome's genes and exon blocks: `og map genes --chrom` |
| `/map/genes/names.json` | every gene's name and place, for the search: `og map genes --names` |

Each file is gzip-wrapped JSON. The map fetches its document once, a chromosome's genes the first time a view of it reaches the region level, and the names the first time someone types in the search box.

### The document

The document holds columns, not rows: for each chromosome, one array per field, rows sorted by start. Starts are stored as the gap from the previous row's start (`ds`) and ends as lengths (`l`), so the page rebuilds the starts, the ends and their running maximum in one pass. Every span is 0-based with an exclusive end, as BED; a single position is 1-based, as VCF.

```jsonc
{
  "v": 1, "build": "GRCh38", "kind": "person", "run": "DEMO-001",
  "people": [{"id": "DEMO-001", "label": "DEMO-001", "sex": "male", "has_y": true}],
  "ref": {"chroms": [{"n": "chr1", "len": 248956422, "cen": [121700000, 125100000], "masked": [[121700000, 143200000]]}],
          "bands": {"chr1": {"ds": [], "l": [], "name": [], "stain": []}}, "gencode": "v50"},
  "layers": {"findings": {"state": "ok", "count": 9}, "roh": {"state": "ok", "count": 2}},
  "per": {"DEMO-001": {"findings": {"id": [], "level": [], "chrom": [], "s": [], "e": [], "gene": [], "veil": []},
                       "roh": {"chr1": {"ds": [], "l": [], "long": []}}}},
  "dict": {"veil": []},
  "text": {"genome": [], "chrom": {}, "findings": {}}
}
```

A gene file holds `g` (each gene's start, length, strand, name and type) and `x` (the merged exon blocks), in the same gap and length columns. `rust/crates/og-map/README.md` gives every field.

## Speed

The map draws on one canvas. The genome rows are drawn once and kept, each lane finds the rows in view by a binary search, and text widths are measured once. `tests/test_live_demo.py` checks the budget in a headless Chrome on the demo with a genome-scale document (about 4,000 runs of homozygosity and a chromosome with 5,000 genes): the first paint within 150 ms of the document arriving (parsed within 30 ms), and over a scripted pan of 300 steps through the chromosome, region and gene levels, the work of each frame (the input handler and the drawing) at most 8 ms at the median and 16 ms at the 95th percentile. The same check run with a 20 ms busy loop in the drawing must fail.
