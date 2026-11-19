# Bring Your Results into OwnHealth

## What This Does

[OwnHealth](https://github.com/GeiserX/OwnHealth) is a self-hosted health record for a household. It keeps lab results, documents, wearables and, for a genome, a short list of findings next to the rest of a person's record. OwnHealth runs no genome analysis. It reads one file that OwnGenome writes:

```text
${GENOME_DIR}/${SAMPLE}/facts/${SAMPLE}_genome_facts.json
```

`run-all.sh` writes it at the end of every run, after the reports. If step 24 (the HTML report) fails, it skips the file, because `summary.json` may then be an earlier run's. If `bin/og facts` fails, there is no file. `run-all.sh` removes an earlier run's file before the pipeline starts, so a failed or interrupted run never leaves a stale one to import by mistake. A run that stops at its own checks, for example on an unknown step in `TOOLS`, changes nothing: the earlier file stays, and it still matches the earlier run's results. The run's output says why the file was not written. A plain `nextflow run main.nf` does not write it yet: run the `bin/og facts` command under [Write the file by hand](#write-the-file-by-hand) on the folder the run published. The format is `ownhealth.genome-facts/v1`. OwnHealth owns that format and checks every file against it when you import. A file it refuses is not stored at all.

## Bring the file into OwnHealth

1. In OwnHealth, open **Import & review** and pick the person in the household list.
2. Click **Upload**, then **Choose files**, and select these files from the sample folder:
   - `facts/${SAMPLE}_genome_facts.json`
   - `vcf/${SAMPLE}.vcf.gz` and `vcf/${SAMPLE}.vcf.gz.tbi`
   - optional: `vcf/${SAMPLE}.g.vcf.gz` and `vcf/${SAMPLE}.g.vcf.gz.tbi`

   **Choose an OwnGenome folder** does the same in one step: pick a folder that holds only these files, and OwnHealth sends every `.json`, `.vcf.gz` and `.tbi` in it.
3. Wait until the facts file shows **Imported** and the VCF files show **Attached**.
4. Open **Genome**.

The facts file names each VCF by its SHA-256, and OwnHealth keeps a copy that matches. With the VCF, OwnHealth can answer a variant lookup at any position. With the gVCF too, it can also tell "you match the reference here" from "this position was not read". Without them the import still works, but a lookup can only answer the positions the findings cite.

When you run the pipeline again on the same sample, import the new file. It replaces the older one for that sample. OwnHealth refuses a facts file that is already imported for another person.

## What the file holds

| Kind | From | What OwnHealth shows |
|---|---|---|
| PGx | PharmCAT ([step 7](07-pharmacogenomics.md)) and the CPIC lookup ([step 27](27-cpic-lookup.md)), CYP2D6 from the consensus of [step 36](36-pgx-consensus.md) | diplotype, phenotype, the drugs with a CPIC or DPWG guideline |
| HLA | the drug-related alleles of [step 8](08-hla-typing.md) | present or absent |
| Monogenic, ACMG SF, carrier | the Pathogenic and Likely pathogenic hits of the ClinVar screen ([step 6](06-clinvar-screen.md)); a hit is filed under the gene the ACMG SF list reports it for. A hypomorph or pseudodeficiency allele of the carrier list is not a carrier finding on its own, so it is left out and its gene is partly assessed. A gene whose X-linked condition is dominant (MECP2, for example) is never a carrier finding: one change can affect a woman | the variant, its ClinVar review stars, and "a clinical test confirms this before any decision" |
| PRS | the polygenic scores of [step 25](25-prs.md) | the percentile, the PGS Catalog id, the share of the score's variants found |
| Risk markers | APOE from the traits panel ([step 42](42-traits.md)), Lp(a) from [step 46](46-lpa.md) | the genotype or the estimate |

Every finding carries where it came from (the tool and the database version), an evidence level and a confidence level:

- **Evidence:** E1 is a CPIC guideline gene or a ClinVar hit with 2 or more review stars. E2 is a ClinVar hit with 1 star or a replicated association. E3 is a single study or a polygenic score.
- **Confidence:** high when the call passed its quality checks, medium or low when it did not, with the reason. A PGx call drops to medium when a site that defines it has a low genotype quality. CYP2D6 is high only when two callers agree. A run that started from a VCF, with no reads, caps a ClinVar hit at medium and CYP2D6 at low.
- **Call state:** a variant finding, a PGx call and an HLA result say whether they rest on the reference (`ref`) or on a variant (`alt`). A PGx call with an allele PharmCAT could not call, such as `*1/Unknown`, is `no_call`. A polygenic score has no call state. A PGx call is `ref` when both alleles are the gene's reference allele as PharmCAT names it (`*1`, `Reference`, CYP2C19 `*38`, CFTR `No CPIC variants found`). ABCG2, IFNL3 and VKORC1 also cite their one site with its genotype.

The file also holds a coverage map. It lists every gene of the ACMG SF list, the CPIC genes and the carrier panel, and any other gene with a PGx or HLA call (CYP3A4, for example). Each gene is marked as assessed, partly assessed or not assessed. A gene that is not fully assessed has a reason, for example low depth or a near copy of the region elsewhere in the genome. OwnHealth shows the map under the findings, so the reader can tell "no finding" from "not looked at".

The ClinVar screen counts a gene as assessed only when the run shows its reads were there: gene callability ([step 16b](16b-mosdepth.md)) ran and the genome's mean depth reached its floor of 10 reads. Without that, for example on a run from a VCF, its genes are not assessed. An HLA gene is assessed when T1K typed it, and not assessed when T1K did not; step 16b and the mean depth do not apply to it. A call that the file leaves out (see below) marks its gene as partly assessed, with a note that names the call. A gene with a PGx or HLA call is never marked not assessed: when the ClinVar screen cannot show the gene was read, the gene is partly assessed, because the call stands on the sites that define it. A PGx call where PharmCAT could not call one allele also marks its gene as partly assessed. So does a PGx call that step 27 marks ambiguous (possible diplotypes with different phenotypes): it is left out, because it is no single result. A PGx call can stand on its sites while step 16b or the carrier screen cannot show the rest of the gene was read. The gene is then partly assessed, and its row keeps their reason and callable share, for example "paralog" and 41%.

## What OwnHealth does not show

- **The rest of the reports.** Structural variants, repeat expansions, copy numbers, mitochondria, runs of homozygosity, ancestry, haplogroups and the plain traits stay in the [HTML report](24-html-report.md). OwnHealth's format has no place for them yet.
- **Sensitive results, by default.** APOE, the ACMG SF cancer and cardiac genes and predictive neurological repeats are tagged as sensitive or adult onset. OwnHealth shows them as a locked row until the person chooses to reveal each one. A separate setting decides whether an agent may see them.
- **Adult-onset results for children.** For a person under 18, OwnHealth stores some results but does not show them. These are adult-onset and sensitive results, carrier status, and scores for adult-onset conditions. At 18 the person decides.
- **Words outside OwnHealth's terms.** A PGx phenotype or diplotype that OwnHealth's word list does not hold (for example "Uncertain Susceptibility") is left out of the finding; the call and its call state stay. PharmCAT's placeholders for "no phenotype" (`n/a`, `No Result`, `Indeterminate`) are left out the same way. The file keeps each call's `value.call_state`, but OwnHealth today takes the call state only from the sites a call cites or from a `*1` diplotype. Until it reads `value.call_state`, it can show a reference call such as CACNA1S `Reference/Reference` as `alt` and actionable.
- **Findings the file cannot place.** `og facts` leaves out a variant whose reference base does not match the reference FASTA, a variant whose reference bases or neighbours fall in a masked stretch of the reference, and an HLA allele T1K could not type. The coverage map marks the gene as partly assessed, and the log of the step, `logs/ownhealth-facts.log`, names each one.

The findings are information for a personal record, not a diagnosis. A ClinVar hit is a lead to confirm with a clinical test, never a result to act on alone.

## Write the file by hand

`run-all.sh` runs this for you. To write the file again after a single step, first run step 24 again, because `og facts` reads the results from `summary.json` and only step 24 (or `generate-report.sh`) writes it fresh. Remove the earlier file first, so a failed step 24 leaves no file to import:

```bash
rm -f "${GENOME_DIR}/${SAMPLE}/facts/${SAMPLE}_genome_facts.json"
./scripts/24-html-report.sh "${SAMPLE}" && ./scripts/ownhealth-facts.sh "${SAMPLE}"
```

It needs `summary.json` (step 24 or `generate-report.sh` writes it), `og` (`./scripts/setup.sh --og` builds it, see [The Rust core](rust.md)) and the GRCh38 reference FASTA the run used, with its `.fai`. The script runs `og` in Docker, in the image `PYTHON_IMAGE` of `versions.env`. It takes a few seconds for a whole genome: it reads `summary.json`, the PharmCAT report and the ClinVar hits, a few reference bases per variant, and hashes the VCF and gVCF.

The variant identifiers in the file are [GA4GH VRS](https://vrs.ga4gh.org/) 2.0 identifiers worked out on that reference. `og facts` refuses a reference that names none of GRCh38's chromosomes (`chr1` to `chrM`) or has one of a length other than GRCh38's. It cannot check every base, so give it the reference the run used.

The same command without the script, on the host. Give `--reference` the FASTA the run used: the default below, or the file `REF_FASTA` named for the run:

```bash
bin/og facts --sample "${SAMPLE}" --sample-dir "${GENOME_DIR}/${SAMPLE}" \
  --reference "${REF_FASTA:-${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta}"
```

Exit code 0 means the file was written, 1 a missing or unreadable input, 2 a wrong argument. The file is readable by its owner only.

## Files

| File | What it is |
|---|---|
| `assets/ownhealth/genome-facts.v1.json` | OwnHealth's schema for the file, copied unchanged |
| `assets/ownhealth/sensitive_genes.tsv` | OwnHealth's list of sensitive genes and their positions, copied unchanged |
| `assets/ownhealth/VENDORED.tsv` | the OwnHealth commit and SHA-256 of each copied file. A test fails when a copy changes |
| `assets/ownhealth/cpic_genes.tsv` | the genes with a CPIC guideline, for the coverage map |
| `assets/ownhealth/grch38_refget.tsv` | the GA4GH refget identifier of each GRCh38 chromosome, for the VRS identifiers. `scripts/ci/build-refget-table.sh` rebuilds it |
