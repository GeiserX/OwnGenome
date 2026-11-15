# og-map

The genome map: one document per person that the live page draws, built on request from the run folder and the repository's reference assets. This crate holds the rules. The page draws what the document says.

This crate today: the reference bands, the chromosome lengths, centromeres and masked spans, the findings layer with the veil tags, the runs-of-homozygosity layer, and the gene files. Structural variants, mosaic calls and ClinVar hits come in later layers, each as a new key; a reader ignores keys it does not know.

## The verbs

```bash
og map build --root GENOME_DIR --run RUN --assets ASSETS --cache-dir DIR < facts.json > map.json.gz
```

- `--root`: the folder with one folder per run (GENOME_DIR). `RUN` is a run name as the live site accepts it (a letter or digit, then letters, digits, `.`, `_` or `-`).
- `--assets`: the repository's `assets/`. The bands come from `assets/reference/cytoBand.hg38.txt` (`assets/reference/README.md` says where it comes from).
- `--cache-dir`: where built documents are kept, as `DIR/RUN/<stamp>.json.gz`. A new stamp removes the run's older ones.
- stdin: `{"results": <the results body the page gets>}`, or nothing. The findings come from `results.findings`.
- The person's sex is the one step 16 infers in `RUN/RUN_sex_check.tsv` when it says male or female, else `results.person.inferred_sex`, else `declared_sex`, else unknown.
- Gene-only findings are placed with the highest GENCODE release installed under `GENOME_DIR/reference` (the genes file `og map genes` reads). Without one they stay off the map, and `ref.gencode` is `null`.
- stdout: the document.
- `--predicate-table`: instead of the document, what the map draws one by one at each level, for every row of the document and every gene of the installed GENCODE release: a tab-separated header `layer chrom s e id genome chromosome region gene reads`, then one line per row with `1` where it draws one by one and `0` where it does not. The rule is `og_map::predicate::draws_at`, the same as `drawsAt()` in `assets/live/map.js`: a finding at every level, a run of homozygosity of 5 Mb or more at every level and a shorter one from the chromosome level, a gene from the region level. A closed finding is listed too: the veil is the page's rule. The cache is neither read nor written. The live demo's browser check compares the two row by row.

Exit codes: 0 the document is on stdout; 2 usage (bad option, bad run name, stdin that is not a JSON object); 3 no map for this run (no folder for it under the root), with the reason on stderr; 1 an internal error. A layer that cannot be read is not an exit code: it is that layer's `unreadable` state.

**The stamp** is a 64-bit FNV-1a of the format version, the og tree hash, the stdin bytes, and the path, size and modification time of every input, present or absent: the band file, the sex file, the GENCODE genes file (or, while none is installed, the `reference` folder), and the ROH texts. It detects change; it is not a defence against someone who wants a collision. A cache hit returns the stored bytes without reading any input.

`og map genes` writes the gene files the map draws, apart from any run:

```bash
og map genes --root GENOME_DIR --cache-dir DIR --chrom chr7 > genes.chr7.json.gz
og map genes --root GENOME_DIR --cache-dir DIR --names     > genes.names.json.gz
```

- Inputs: `GENOME_DIR/reference/gencode.v<N>.basic.genes.gtf` (the gene lines) and `gencode.v<N>.basic.exons.bed` (every exon of a basic transcript), which setup installs. When several releases are installed, the highest `N` is read.
- `genes.<chrom>.json.gz`: `{"v":1,"gencode":"v50","chrom":"chr7","g":{"ds":[],"l":[],"strand":[],"name":[],"type":[]},"x":{"ds":[],"l":[]}}`. `g` holds every gene on that chromosome, `x` the exon blocks: overlapping or touching exons merged, so a block is "an exon of some basic transcript here", not one transcript. Rows are sorted by start, as gaps and lengths. Only chr1-22, X and Y are kept.
- `genes.names.json.gz`: `{"v":1,"gencode":"v50","name":[],"chrom":[],"s":[],"e":[]}`, every gene sorted by name, with absolute spans.
- GTF positions are 1-based and inclusive, so `s = start - 1` and `e = end`; the BED is already 0-based.
- The first call builds all 25 files and caches them under `DIR/genes/<stamp>/` (the stamp covers both inputs' path, size and modification time, the format version and the og tree hash). Later calls stream the cached file. A new set removes the older ones.
- Exit codes: 2 for usage (also a chromosome not on the map); 3 when the GENCODE files are missing or cannot be read, the reason starting "Gene names need step 08's GENCODE file".

## The document, version 1

One gzip member with mtime 0 and no file name, holding JSON written as Python's `json.dumps(doc, separators=(",", ":"))` with the keys in the order below. The same input always gives the same bytes.

- Every span is 0-based with an exclusive end, as BED. A single position `pos` is 1-based, as VCF.
- Columnar layers: per chromosome, an object of equal-length arrays, rows sorted by start. Starts are stored as `ds`, the gap from the previous row's start (the first row's value is absolute), and ends as `l`, the length (`e - s`). A reader rebuilds `s`, `e` and the running maximum of `e` in one pass.

```jsonc
{
  "v": 1, "build": "GRCh38", "kind": "person",
  "run": "DEMO-001",
  "stamp": "<16 hex digits>",
  "made_by": "<og tree hash>",
  "people": [{"id": "DEMO-001", "label": "DEMO-001", "sex": "male|female|unknown", "has_y": true}],  // has_y null when the sex is unknown
  "ref": {
    "chroms": [{"n": "chr1", "len": 248956422, "cen": [121700000, 125100000], "masked": [[121700000, 143200000]]}],
    "bands": {"chr1": {"ds": [], "l": [], "name": [], "stain": []}},
    "gencode": "v50"                   // null when no GENCODE release is installed
  },
  "layers": {"findings": L, "roh": L},
  "per": {
    "DEMO-001": {
      "findings": {"id": [], "level": [], "gate": [], "chrom": [], "pos": [], "s": [], "e": [], "gene": [], "title": [],
                   "hgvs": [], "zyg": [], "copy": [], "steps": [], "href": [], "reads": [], "phrase": [], "veil": []},
      "roh": {"chr1": {"ds": [], "l": [], "markers": [], "q": [], "long": [], "art": []}},
      "roh_bins": {"size": 1000000, "chr1": []}
    }
  },
  "dict": {"art": ["", "centromere", "one X", "one Y"], "veil": ["<finding id>"]},
  "text": {"genome": [{"t": "2 runs", "b": true}, {"t": " of homozygosity of 5 Mb or more, 12.5 Mb in all.", "b": false}],
           "chrom": {"chr1": [{"t": "…", "b": false}]},
           "findings": {"none": [], "one": [], "many": [{"t": "{n} findings", "b": true}, {"t": " have a place on the map.", "b": false}]}}
}
```

- `ref.chroms`: chr1-22, X and Y. `len` is the chromosome's last band end, `cen` the span of its `acen` bands, `masked` its `acen`, `gvar` and `stalk` bands with touching neighbours merged.
- `L`, a layer entry: `state` (`ok`, `missing` or `unreadable`), `step`, `step_name`, `source` (the file read or looked for, relative to GENOME_DIR), `note` (why it is missing, or `file:line: reason` when unreadable), `count` (rows on the map), `totals` (`all` rows read and `kept` rows on the map, when the file was read), `caveat`, and `empty` (the sentence for a read file with no rows). A missing input is `missing`, an input that cannot be read is `unreadable`, and a read input with no rows is `ok` with `count` 0: three answers, never a blank lane.
- `findings`: one flat column per field, rows sorted by chromosome then start. A finding's `viewer` (`chrom:pos` or `chrom:start-end`, 1-based, `chr` optional) places it: a position gives `pos` and `[pos - 1, pos)`, a span gives `[start - 1, end)` and `pos` null. Without a usable `viewer`, a finding that names a gene sits on that gene's span from the GENCODE genes file (the first in map order when a name is on several chromosomes). A finding with neither is not a row. `gate` is set for a finding the page keeps behind its own warning: the traits card holds APOE that way, so the map adds the finding `apoe` on the APOE gene, with `gate` `apoe`; the hidden text, which names the genotype, never enters the document. `zyg` (`het`, `hom`, `hemi` or `unknown`) comes from the finding's own `Genotype` evidence; no VCF is read. `hgvs` is its `Variant` evidence. `copy` is null until the family map. `href` is the finding's page relative to the person's page, `reads` the place as the reads viewer takes it (1-based), `phrase` the sentence the map shows. `veil` is the index in `dict.veil` of the first closable finding (in `dict.veil` order) whose span the row overlaps: its own for a closable finding with no earlier overlap, else -1. So an open finding over a closable one's span is tagged too. The layer's `count` and `totals.kept` are the rows whose `veil` is -1; `totals.all` adds the findings that are not closable and have no place, and the `note` counts those ("N findings with no place on the map stay on the page").
- **The veil.** The document is the same for every reader. Whether a page is closed is the page's own rule; the document says only what a closed result could reveal. Closable findings are those of level `act` and every gated one (APOE closes on every page, an adult's too). Each closable finding with a place veils its span, and `dict.veil` lists their ids: gated ones first, then the act ones in the page's order. A row of a layer that can be a finding's own source (structural variants, mosaic calls and ClinVar hits, as their layers arrive), and every findings row, that overlaps a veiled span carries the index of the first such finding in `dict.veil` in a `veil` column, else -1. A row holds one tag only. Gated findings come first, so a row over both APOE and an act finding carries APOE, which is closed on every page; a row over two act findings carries the first, and both close under the same rule. The page decides its own findings with its own veil, and the tag binds only the rows that are not one of the page's findings (APOE, or a row the page's results do not have yet): the page hides such a row while its result is closed. Counts, totals and density bins cover untagged rows only, so no number counts a closed row. Runs of homozygosity and bands are never a finding's source and carry no tag.
- `roh`: the `RG` rows of step 11's `bcftools roh` text (`{run}/vcf/{run}_roh.txt`, else `{run}/roh/{run}_roh.txt`). bcftools writes 1-based inclusive positions, so `s = start - 1` and `e = end`. `ST` lines are skipped; rows on chrM and other contigs count in `totals.all` only. `markers` and `q` are bcftools' columns 7 and 8. `art` indexes `dict.art`: `one X` on a male's chrX, `one Y` on a male's chrY, else `centromere` when more than half the run lies in masked spans. An unknown sex sets neither flag. `long` is 1 for a run of 5 Mb or more with no `art`.
- `roh_bins`: bases covered by any run per 1 Mb bin, for each chromosome with a run.
- `text`: the sentences the page shows, in pieces; `b` marks bold words. `text.findings` is the count of findings on the map in three forms (`none`, `one`, `many`): the page counts the rows it draws: each of its own findings that its veil shows at this moment, whatever the row's tag, and each other row whose tag is -1 or names a result that is open; it picks the form and puts the number in place of `{n}`.

## Tests

```bash
cargo test -p og-map
cargo run --release -p og-map --example genome_scale -- OUT_DIR [SEED]   # a seeded, invented 200 MB ROH text with about 4,000 runs
```

`tests/veil.rs` and `tests/findings.rs` read a synthetic person under `tests/data/facts` (invented findings and rows, a `DEMO-005` id, no genotype) with the public GENCODE slice under `tests/fixtures/map/gencode-slice`.

`tests/contract.rs` fails when a path og-map reads is not a `file` row of `assets/output_contract.tsv`, so a renamed output breaks the build loudly instead of turning a layer to `missing`.
