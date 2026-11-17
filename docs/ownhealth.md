# Bring Your Results into OwnHealth

## What This Does

[OwnHealth](https://github.com/GeiserX/OwnHealth) is a self-hosted health record for a household. It keeps lab results, documents, wearables and, for a genome, a short list of findings next to the rest of a person's record. OwnHealth runs no genome analysis. It reads one file that OwnGenome writes:

```text
${GENOME_DIR}/${SAMPLE}/facts/${SAMPLE}_genome_facts.json
```

`run-all.sh` writes it at the end of every run, after the reports. The format is `ownhealth.genome-facts/v1`. OwnHealth owns that format and checks every file against it when you import. A file it refuses is not stored at all.

## Bring the file into OwnHealth

1. In OwnHealth, open **Import & review** and pick the person in the household list.
2. Click **Choose files** and select these files from the sample folder:
   - `facts/${SAMPLE}_genome_facts.json`
   - `vcf/${SAMPLE}.vcf.gz` and `vcf/${SAMPLE}.vcf.gz.tbi`
   - optional: `vcf/${SAMPLE}.g.vcf.gz` and `vcf/${SAMPLE}.g.vcf.gz.tbi`
3. Wait until the facts file shows **Imported** and the VCF files show **Attached**.
4. Open **Genome**.

The facts file names each VCF by its SHA-256, and OwnHealth keeps a copy that matches. With the VCF, OwnHealth can answer a variant lookup at any position. With the gVCF too, it can also tell "you match the reference here" from "this position was not read". Without them the import still works, but a lookup can only answer the positions the findings cite.

When you run the pipeline again on the same sample, import the new file. It replaces the older one for that sample. OwnHealth refuses a facts file that is already imported for another person.

## What the file holds

| Kind | From | What OwnHealth shows |
|---|---|---|
| PGx | PharmCAT ([step 7](07-pharmacogenomics.md)) and the CPIC lookup ([step 27](27-cpic-lookup.md)), CYP2D6 from the consensus of [step 36](36-pgx-consensus.md) | diplotype, phenotype, the drugs with a CPIC or DPWG guideline |
| HLA | the drug-related alleles of [step 8](08-hla-typing.md) | present or absent |
| Monogenic, ACMG SF, carrier | the Pathogenic and Likely pathogenic hits of the ClinVar screen ([step 6](06-clinvar-screen.md)) | the variant, its ClinVar review stars, and "a clinical test confirms this before any decision" |
| PRS | the polygenic scores of [step 25](25-prs.md) | the percentile, the PGS Catalog id, the share of the score's variants found |
| Risk markers | APOE from the traits panel ([step 42](42-traits.md)), Lp(a) from [step 46](46-lpa.md) | the genotype or the estimate |

Every finding carries where it came from (the tool and the database version), an evidence level and a confidence level:

- **Evidence:** E1 is a CPIC guideline gene or a ClinVar hit with 2 or more review stars. E2 is a ClinVar hit with 1 star or a replicated association. E3 is a single study or a polygenic score.
- **Confidence:** high when the call passed its quality checks, medium or low when it did not, with the reason. A PGx call drops to medium when a site that defines it has a low genotype quality. CYP2D6 is high only when two callers agree.

The file also holds a coverage map. It lists every gene of the ACMG SF list, the CPIC genes and the carrier panel. Each gene is marked as assessed, partly assessed or not assessed. A gene that is not fully assessed has a reason, for example low depth or a near copy of the region elsewhere in the genome. OwnHealth shows the map under the findings, so the reader can tell "no finding" from "not looked at".

## What OwnHealth does not show

- **The rest of the reports.** Structural variants, repeat expansions, copy numbers, mitochondria, runs of homozygosity, ancestry, haplogroups and the plain traits stay in the [HTML report](24-html-report.md). OwnHealth's format has no place for them yet.
- **Sensitive results, by default.** APOE, the ACMG SF cancer and cardiac genes and predictive neurological repeats are tagged as sensitive or adult onset. OwnHealth shows them as a locked row until the person chooses to reveal each one. A separate setting decides whether an agent may see them.
- **Adult-onset results for children.** For a person under 18, OwnHealth stores some results but does not show them. These are adult-onset and sensitive results, carrier status, and scores for adult-onset conditions. At 18 the person decides.
- **Findings that do not fit OwnHealth's terms.** `og facts` leaves out a finding that it cannot describe with the words OwnHealth accepts, for example a PGx phenotype outside OwnHealth's list. It also leaves out a variant whose reference base does not match the reference FASTA. The log of the step, `logs/ownhealth-facts.log`, names each one.

The findings are information for a personal record, not a diagnosis. A ClinVar hit is a lead to confirm with a clinical test, never a result to act on alone.

## Write the file by hand

`run-all.sh` runs this for you. To write the file again after a single step:

```bash
./scripts/ownhealth-facts.sh ${SAMPLE}
```

It needs `summary.json` (step 24 or `generate-report.sh` writes it), `og` (`./scripts/setup.sh --og` builds it, see [The Rust core](rust.md)) and the reference FASTA the run used. The variant identifiers in the file are [GA4GH VRS](https://vrs.ga4gh.org/) 2.0 identifiers worked out on that reference, so a reference that is not GRCh38 is refused.

The same command without the script, on the host:

```bash
bin/og facts --sample ${SAMPLE} --sample-dir ${GENOME_DIR}/${SAMPLE} \
  --reference ${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta
```

Exit code 0 means the file was written, 1 a missing or unreadable input, 2 a wrong argument.

## Files

| File | What it is |
|---|---|
| `assets/ownhealth/genome-facts.v1.json` | OwnHealth's schema for the file, copied unchanged |
| `assets/ownhealth/sensitive_genes.tsv` | OwnHealth's list of sensitive genes and their positions, copied unchanged |
| `assets/ownhealth/VENDORED.tsv` | the OwnHealth commit and SHA-256 of each copied file. A test fails when a copy changes |
| `assets/ownhealth/cpic_genes.tsv` | the genes with a CPIC guideline, for the coverage map |
| `assets/ownhealth/grch38_refget.tsv` | the GA4GH refget identifier of each GRCh38 chromosome, for the VRS identifiers. `scripts/ci/build-refget-table.sh` rebuilds it |
