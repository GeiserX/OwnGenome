# og-map

The genome map: one document per person that the live page draws, built on request from the run folder and the repository's reference assets. This crate holds the rules. The page draws what the document says.

Phase 1 core (this crate today): the reference bands, the chromosome lengths, centromeres and masked spans, and the runs-of-homozygosity layer. The findings, structural variants, mosaic calls, ClinVar hits and genes come in later layers, each as a new key; a reader ignores keys it does not know.

## The verb

```bash
og map build --root GENOME_DIR --run RUN --assets ASSETS --cache-dir DIR < facts.json > map.json.gz
```

- `--root`: the folder with one folder per run (GENOME_DIR). `RUN` is a run name as the live site accepts it (a letter or digit, then letters, digits, `.`, `_` or `-`).
- `--assets`: the repository's `assets/`. The bands come from `assets/reference/cytoBand.hg38.txt` (`assets/reference/README.md` says where it comes from).
- `--cache-dir`: where built documents are kept, as `DIR/RUN/<stamp>.json.gz`. A new stamp removes the run's older ones.
- stdin: `{"results": <the results body the page gets>}`, or nothing. The person's sex is `results.person.inferred_sex` when it is male or female, else `declared_sex`, else unknown.
- stdout: the document.

Exit codes: 0 the document is on stdout; 2 usage (bad option, bad run name, stdin that is not a JSON object); 3 no map for this run (no folder for it under the root), with the reason on stderr; 1 an internal error. A layer that cannot be read is not an exit code: it is that layer's `unreadable` state.

**The stamp** is a 64-bit FNV-1a of the format version, the og tree hash, the stdin bytes, and the path, size and modification time of every input, present or absent. It detects change; it is not a defence against someone who wants a collision. A cache hit returns the stored bytes without reading any input.

## The document, version 1

One gzip member with mtime 0 and no file name, holding JSON written as Python's `json.dumps(doc, separators=(",", ":"))` with the keys in the order below. The same input always gives the same bytes.

- Every span is 0-based with an exclusive end, as BED.
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
    "gencode": null
  },
  "layers": {"roh": L},
  "per": {
    "DEMO-001": {
      "roh": {"chr1": {"ds": [], "l": [], "markers": [], "q": [], "long": [], "art": []}},
      "roh_bins": {"size": 1000000, "chr1": []}
    }
  },
  "dict": {"art": ["", "centromere", "one X"]},
  "text": {"genome": [{"t": "2 runs", "b": true}, {"t": " of homozygosity of 5 Mb or more, 12.5 Mb in all.", "b": false}],
           "chrom": {"chr1": [{"t": "…", "b": false}]}}
}
```

- `ref.chroms`: chr1-22, X and Y. `len` is the chromosome's last band end, `cen` the span of its `acen` bands, `masked` its `acen`, `gvar` and `stalk` bands with touching neighbours merged.
- `L`, a layer entry: `state` (`ok`, `missing` or `unreadable`), `step`, `step_name`, `source` (the file read or looked for, relative to GENOME_DIR), `note` (why it is missing, or `file:line: reason` when unreadable), `count` (rows on the map), `totals` (`all` rows read and `kept` rows on the map, when the file was read), `caveat`, and `empty` (the sentence for a read file with no rows). A missing input is `missing`, an input that cannot be read is `unreadable`, and a read input with no rows is `ok` with `count` 0: three answers, never a blank lane.
- `roh`: the `RG` rows of step 11's `bcftools roh` text (`{run}/vcf/{run}_roh.txt`, else `{run}/roh/{run}_roh.txt`). bcftools writes 1-based inclusive positions, so `s = start - 1` and `e = end`. `ST` lines are skipped; rows on chrM and other contigs count in `totals.all` only. `markers` and `q` are bcftools' columns 7 and 8. `art` indexes `dict.art`: `one X` on a male's chrX, else `centromere` when more than half the run lies in masked spans. An unknown sex sets no `one X` flag. `long` is 1 for a run of 5 Mb or more with no `art`.
- `roh_bins`: bases covered by any run per 1 Mb bin, for each chromosome with a run.
- `text`: the sentences the page shows, in pieces; `b` marks bold words.

## Tests

```bash
cargo test -p og-map
cargo run --release -p og-map --example genome_scale -- OUT_DIR [SEED]   # a seeded, invented 200 MB ROH text with about 4,000 runs
```

`tests/contract.rs` fails when a path og-map reads is not a `file` row of `assets/output_contract.tsv`, so a renamed output breaks the build loudly instead of turning a layer to `missing`.
