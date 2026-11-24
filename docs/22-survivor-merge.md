# Step 22: Structural Variant Consensus Merge

## What This Does

Keeps the structural variants (SVs) that two or more callers agree on. SURVIVOR merge pairs two calls when both of their breakpoints lie within 1,000 bp of each other, their SV type and strands agree, and the event is at least 50 bp long. Then `og sv consensus-check` keeps a merged record only when its calls are one event: the same type, the same strands, and an overlap of at least half of the longer call (see [One event per record](#one-event-per-record)). Each consensus record says how many callers support it and which. That is caller agreement, not validation (see Why).

A caller votes only with the calls it makes (see [Which records vote](#which-records-vote)): a record the caller genotyped as reference, or filtered, is not a second caller.

## Why

Individual SV callers each have distinct biases and false-positive profiles:

- **Manta**: fast, sensitive for smaller SVs (paired-end and split reads)
- **Delly**: strongest for inversions and balanced translocations (paired-end, split reads and depth)
- **CNVpytor**: best for large CNVs (read depth only)

An SV that two callers report is less likely to be one caller's artefact. Agreement is not independent proof, though. Manta and Delly both read split and paired reads from the same alignment, so a misaligned read cluster both of them see can reach the consensus. Only CNVpytor reads a different signal, depth, and its breakpoints are coarse (see Limitations).

SURVIVOR's 1,000 bp rule alone is not enough. The pinned SURVIVOR 1.0.7 merges two deletions of 100 bp that start 900 bp apart and share no base, and on some inputs it merges a deletion and a duplication at the same place into one record, although its type and strand flags are set. Each such record read "2 callers" for an event no two callers called. `og sv consensus-check` reads the calls inside each record and splits it by event, so a record means one event again. The Rust tests of `og-steps` (`rust/crates/og-steps/src/sv/tests.rs`) run SURVIVOR's own records of these cases, and the e2e case `tests/e2e/sv-merge-consensus-check.sh` runs them through the step.

The step used to bin calls by chromosome, `int(POS/1000)` and SV type. That split one deletion called at positions 999 and 1001 over two bins, and counted two deletions of very different size that start in one bin as agreement. SURVIVOR compares both breakpoints, so neither happens: the e2e case `tests/e2e/sv-mito-telomere-steps-1-sv-merge.sh` runs both situations on the synthetic three-caller set in `tests/fixtures/sv/`. The e2e case `tests/e2e/sv-merge-votes.sh` runs the records that must not vote on the same set.

## Tool

- **SURVIVOR** 1.0.7 (Jeffares et al., Nat Commun 2017), `SURVIVOR merge`
- **og** (`og sv survivor-prep`) for the records each caller votes with, and (`og sv consensus-check`) for the check that each merged record is one event
- **bcftools** for the sorted, indexed output

Not chosen: `truvari collapse`, which also compares sequence similarity. It would need a comparison run against SURVIVOR on real calls first.

## Docker Image

- `SURVIVOR_IMAGE` (the merge), `BCFTOOLS_IMAGE` (runs `og` for the records that vote and for the check of each merged record; sort, index). In Nextflow the check runs in the `SURVIVOR_MERGE` task, right after the merge

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Input

At least two of the following (the script uses every one it finds):

| Caller | Expected path |
|---|---|
| Manta (step 4) | `${GENOME_DIR}/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz` |
| Delly (step 19) | `${GENOME_DIR}/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz` |
| CNVpytor (step 18) | `${GENOME_DIR}/${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.vcf.gz`, or `_cnvs.txt`, which the script turns into that VCF first |
| TIDDIT (script 4a) | `${GENOME_DIR}/${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz` |
| Sniffles2 (script 4c, long reads) | `${GENOME_DIR}/${SAMPLE}/sv_sniffles/${SAMPLE}_sv.vcf.gz` |

GRIDSS (step 4b) is left out of the SV consensus, and the script says so when its VCF exists. GRIDSS reports every event as a pair of breakends (`SVTYPE=BND`), and nothing here turns a pair into one DEL, DUP or INV. SURVIVOR 1.0.7 counts BND as the same type as INV and TRA, so a GRIDSS breakend would pair with another caller's inversion or translocation and count as a vote for an event GRIDSS never typed. Step 4b writes GRIDSS's own high-confidence calls to `sv_gridss/<sample>_gridss.high.vcf.gz`.

The Nextflow pipeline (`survivor_merge` in `--tools`, with at least two of `manta`, `delly`, `cnvpytor`) merges the callers it ran.

## Command

```bash
./scripts/22-survivor-merge.sh your_name
```

## What the Script Does Internally

1. Turns CNVpytor's table into a VCF when step 18 left no VCF, or when the VCF there has no `pytorP1` (an earlier version of this step made it from the table without the e-value). Each record carries the table's e-value, `q0` and `pN` under the INFO names CNVpytor's own VCF uses (`pytorP1`, `pytorQ0`, `pytorPN`)
2. Writes the records each caller votes with (below) as plain VCF under `sv_merged/inputs/`, with one sample column named after the caller (`og sv survivor-prep --caller NAME IN OUT`, which prints how many records it left out and why). SURVIVOR names its output columns after the input samples, so three inputs that all name the sample would give one name three times
3. Runs `SURVIVOR merge inputs/sv_files.txt 1000 2 1 1 0 50`: maximum breakpoint distance 1,000 bp, support from 2 callers, same type, same strands, no size-scaled distance, minimum size 50 bp
4. Checks that SURVIVOR wrote a VCF (it exits 0 when it cannot open an input)
5. Checks each merged record (`og sv consensus-check --callers LIST IN OUT`, below), which prints how many records were one event, how many it split and how many it dropped
6. Sorts, compresses and indexes the result
7. Prints the count, and the count per type and caller combination

### Which records vote

SURVIVOR counts a caller as support for every record of its input. A record goes in only when it is a call the caller makes:

| Record | Votes | Why |
|---|---|---|
| FILTER `PASS` or `.`, genotype with an ALT allele (`0/1`, `1/1`, `1`) | yes | a call |
| genotype `./.`, or no `GT`, or a sites-only VCF | yes | nothing says the call is absent; the VCF made from CNVpytor's table writes `./.` for every call |
| FILTER other than `PASS` or `.` (Delly's `LowQual`) | no | the caller filtered it |
| sample filter `FT` other than `PASS` or `.` (Delly's per-sample `LowQual`) | no | the caller filtered it for this sample |
| genotype with a called allele and no ALT (`0/0` phased or not, haploid `0`, `0/.`) | no | the caller tested the site and found the reference |
| a CNVpytor call with e-value (`pytorP1`) 0.01 or more, `q0` (`pytorQ0`) 0.5 or more, or either value missing | no | CNVpytor's quality rule, the one step 5 applies when it adds CNVpytor's calls ([step 18](18-cnvpytor.md#filtering)) |

Delly genotypes every site it reports, so a PASS Delly record with `0/0` is common. Before this rule, such a record next to a Manta call read as two callers.

### One event per record

`og sv consensus-check --callers LIST IN OUT` reads the same list of prepared caller VCFs as SURVIVOR. It looks up each member by caller, first contig, both numeric breakpoint positions, and the per-member type. It checks the original records for orientation, insertion length, and the mate contig. SURVIVOR summarizes those fields across a caller column or the cluster.

| Calls | Rule |
|---|---|
| All calls | Same type and orientation. A same-contig `BND` matches `INV`. An intercontig `BND` matches `TRA`. |
| Deletions, duplications, inversions | Reciprocal overlap of at least 0.5. The overlap covers at least half of the longer call. |
| Insertions | Starts within 1,000 bp and the shorter call at least half the length of the longer call. |
| Intercontig breakends | Same pair of contigs in the same order and both breakpoints within 1,000 bp. |

The checker selects the largest group in which every pair agrees and each caller votes once. It breaks ties by caller column order and then by the order of the caller's `CO` entries. It removes the selected members and repeats until no group has two callers. A cluster with exactly one verified member for each voting caller and one compatible group passes through unchanged.

A split record takes `POS`, `END`, `REF`, `ALT`, and `ID` from the lower middle member when ordered by start and end. For an intercontig group, the checker uses a retained bracketed allele if the middle member used a symbolic allele and another member supplies a bracketed allele. It takes `QUAL` from the highest retained quality and rebuilds `SUPP`, `SUPP_VEC`, `SVTYPE`, `SVLEN`, `CHR2`, `CIPOS`, `CIEND`, and `STRANDS` from the selected calls. The caller columns name only the retained source IDs. An unnamed source call gets `.` in `FORMAT/ID`.

The checker drops an unmatched or ambiguous member, or a member whose insertion length is missing. The log counts each unresolved vote. It also counts and drops records with malformed `SUPP` or `SUPP_VEC`. Every output record has at least two verified callers.

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_sv_consensus.vcf.gz` | Consensus SVs supported by 2+ callers, each one event ([One event per record](#one-event-per-record)), with `SUPP` and `SUPP_VEC` in INFO. A header line `##og_sv_consensus_check` states the check |
| `${SAMPLE}_sv_consensus.vcf.gz.tbi` | Tabix index |
| `inputs/sv_files.txt`, `inputs/<caller>.vcf` | The records each caller voted with, which SURVIVOR read, in `SUPP_VEC` order |

All output is written to `${GENOME_DIR}/${SAMPLE}/sv_merged/`.

`SUPP_VEC` has one digit per caller, in the order of `sv_files.txt` (manta, delly, cnvpytor, tiddit, sniffles2, for the ones found): `110` is Manta and Delly. Each sample column holds that caller's genotype and its original record ID.

## Runtime

2 to 10 seconds for SURVIVOR's preparation, merge and sort together, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)).

## Interpreting Results

A typical 30X WGS genome produces:

- **Manta**: 7,000-9,000 SVs (see [step 4](04-structural-variants.md))
- **Delly**: 5,000-15,000 SVs
- **CNVpytor**: 3,000-4,000 CNVs, 1,500-2,000 of them with e-value < 0.01 (see [interpreting results](interpreting-results.md#cnvpytor-results-step-18))

An SV two callers report is less likely to be one caller's artefact, but it is not validated: see Why. The count after the merge has not been measured on a real genome with this version of the step.

SV types in the output: **DEL** (deletion), **DUP** (duplication), **INV** (inversion), **INS** (insertion), **TRA** (a translocation, SURVIVOR's name for a breakend pair between two chromosomes).

### Quick inspection

```bash
source versions.env   # from the repository root
# Consensus SVs by type and by which callers support them
docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
  bcftools query -f '%INFO/SVTYPE\t%INFO/SUPP_VEC\n' \
    /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz | sort | uniq -c | sort -rn
```

## Limitations

- Single-caller SVs are discarded even if they are real. If you suspect a specific SV, check the individual caller outputs directly.
- CNVpytor's calls are depth-only and their breakpoints are coarse (its bin size), so a CNVpytor call often lies more than 1 kb from the Manta or Delly breakpoint of the same event and does not pair with it.
- The CNVpytor VCF made from its table carries no genotype (`./.`); its e-value and `q0` come from the table.
- A breakend (BND) record pairs with a breakend, an inversion (INV) or a translocation (TRA) of another caller: SURVIVOR 1.0.7 lets a BND match those two types, and the check keeps that rule. On one chromosome the pair must still share half of the longer span. That includes a lone Manta inversion junction of 1 Mb or more, which step 4 keeps as its two breakends instead of a whole-span inversion.
- The check splits a record and never joins two. Calls SURVIVOR put in separate records stay separate, even when they would pass the check.

## Notes

- Run this step only after completing at least two of: step 4 (Manta), step 19 (Delly), step 18 (CNVpytor).
- Load the consensus VCF into IGV to look at it, but never run AnnotSV on it. Its sample columns are the callers, not you, and its positions are SURVIVOR's merged ones (or, for a record the check split, those of its middle call). AnnotSV (step 5) annotates the callers' own calls, and the reports join `SUPP` and `SUPP_VEC` onto them ([step 5](05-annotsv.md)).

## Links

- [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
- [Manta](https://github.com/Illumina/manta)
- [Delly](https://github.com/dellytools/delly)
- [CNVpytor](https://github.com/abyzovlab/CNVpytor)
