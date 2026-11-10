# Reference data shipped with the repository

## `cytoBand.hg38.txt`

UCSC's chromosome bands for GRCh38 (hg38): one row per band, tab-separated, with the columns chromosome, start (0-based), end (exclusive), band name and Giemsa stain. The genome map (`rust/crates/og-map`) draws the bands from it and derives the chromosome lengths, the centromeres and the masked spans (the `acen`, `gvar` and `stalk` bands, neighbours merged) from it.

| | |
|---|---|
| Source | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz |
| Retrieved | 2026-11-10 |
| sha256 of the download | `e514e48b0a12a516fd5231b44241a20a282f290bdd6131cadcddc6276ee26082`, the pin `scripts/lib/common.sh` checks for the same download |
| Filter | chr1-22, X and Y only, as setup keeps them (the ALT, random, unplaced and chrM rows are dropped) |
| sha256 of this file | `83fe19148485d67da5783aefd9402de58fcd0b0b0d2a476164abeb522522acd1` |
| Rows | 862 |

Rebuild it with:

```bash
curl -sSfO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
sha256sum cytoBand.txt.gz    # e514e48b0a12a516fd5231b44241a20a282f290bdd6131cadcddc6276ee26082
gzip -cd cytoBand.txt.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|X|Y)$/' > cytoBand.hg38.txt
```

### Terms

UCSC's README for the hg38 database folder (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/README.txt) says:

> All the files and tables in this directory are freely usable for any purpose.

UCSC's conditions of use (https://genome.ucsc.edu/conditions.html) say:

> The sequence and annotation data displayed in the Genome Browser are freely available for any use with the following conditions: Usage restrictions for genome sequence data are noted within the species sections on the Credits page. Some annotation tracks contributed by external collaborators contain data that have specific use restrictions.

The chromosome band track carries no use restriction. Credit: the UCSC Genome Browser Group, University of California, Santa Cruz.
