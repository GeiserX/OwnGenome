# Ancestry engine benchmark

Step 26b uses OwnGenome's two-level painter in place of the academic-only RFMix. On simulated genomes with known ancestry, og was 2.2 percentage points behind RFMix in base-pair accuracy for recent three-way admixture (S2) and 12.0 points ahead for older admixture (S4). Improving recent three-way admixture is a follow-up for 1.1.

## Results with the same references

Both painters used the same phased inputs and the same 816 public references after hold-out, from a panel capped at 200 people per continental group. Each row pools five simulated genomes across chromosomes 9, 15, and 21.

| Scenario | og accuracy | RFMix accuracy | og dosage R² | RFMix dosage R² | og seconds | RFMix seconds | og peak MiB | RFMix peak MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S0 | 0.9981 | 0.9927 | 0.9988 | 0.9954 | 14 | 2663 | 140 | 1280 |
| S1 | 0.9614 | 0.9657 | 0.9704 | 0.9738 | 17 | 2850 | 207 | 1287 |
| S2 | 0.8764 | 0.8983 | 0.8829 | 0.9004 | 22 | 2779 | 216 | 1309 |
| S3 | 0.7227 | 0.6267 | 0.7469 | 0.6618 | 15 | 2123 | 182 | 1272 |
| S4 | 0.7817 | 0.6622 | 0.8049 | 0.6899 | 9 | 2880 | 206 | 1280 |

Accuracy is the fraction of bases where both ancestry copies are correct in either order. Dosage R² is the squared correlation of predicted and true copy counts, pooled over bases and groups. Time is the mean wall time per simulated genome for painting. Memory is the highest container peak among the five genomes. On a 4-CPU Linux host, og took 9 to 22 seconds and 140 to 216 MiB. RFMix took 35 to 48 minutes and 1,272 to 1,309 MiB. These measurements cover the three benchmark chromosomes.

| Scenario | Ancestry used to simulate each genome | T |
|---|---|---:|
| S0 | 100% one group, one genome for each of AFR, WEA, CSA, EAS, and AMR | 0 |
| S1 | 80% AFR, 20% WEA | 6 |
| S2 | 50% WEA, 40% AMR, 10% AFR | 12 |
| S3 | 50% WEA, 50% CSA | 20 |
| S4 | 50% EAS, 50% WEA | 50 |

T is the ancestry switch rate used by the simulator, in events per Morgan. It roughly tracks generations since admixture. The groups are African (AFR), West Eurasian (WEA), Central and South Asian (CSA), East Asian (EAS), and Native American (AMR).

| Scenario | Haplotype accuracy (og / RFMix) | Stretches per genome (truth / og / RFMix) | Missed breakpoints (og / RFMix / true) | False breakpoints (og / RFMix) | Fitted T (og) | Fraction mean absolute error (og / ADMIXTURE) |
|---|---:|---:|---:|---:|---:|---:|
| S0 | 0.9990 / 0.9964 | 6 / 6 / 7 | 0 / 0 / 0 | 4 / 5 | 10.3 | 0.001 / 0.001 |
| S1 | 0.9185 / 0.9230 | 23 / 43 / 29 | 19 / 25 / 87 | 102 / 21 | 14.9 | 0.010 / 0.010 |
| S2 | 0.7578 / 0.7843 | 60 / 79 / 58 | 54 / 88 / 273 | 136 / 56 | 18.8 | 0.030 / 0.030 |
| S3 | 0.7258 / 0.7184 | 73 / 56 / 22 | 215 / 294 / 338 | 130 / 37 | 16.1 | 0.008 / 0.008 |
| S4 | 0.7380 / 0.7185 | 189 / 125 / 49 | 309 / 630 / 918 | 112 / 29 | 34.2 | 0.016 / 0.016 |

Haplotype accuracy uses the better of the two phase alignments on each chromosome. Breakpoints count changes in the pair of groups. A true breakpoint is missed when no call is within 1 cM. A called breakpoint is false when no true one is within 1 cM. In S2, og missed fewer true switches but called more false switches. The continental fractions from og and ADMIXTURE `-P` had the same mean absolute error to three decimal places in all five scenarios.

## Simulation and hold-outs

The benchmark starts from a gnomAD HGDP + 1000 Genomes panel installed with `scripts/setup.sh --ancestry-fine-panel`. The earlier hold-out rule chose people whose population labels alone did not justify calling the S0 truth unadmixed. The revised benchmark uses a fixed list of 19 painting references. Each has a projected own-group fraction of at least 0.99 against the full-panel frequencies, with four people per group except for three Karitiana and Surui references for AMR. Restricting AMR to those populations avoids using people whose own-group score is inflated when the panel's AMR frequencies are dominated by Maya references.

The benchmark removes the held-out people from the phasing panel, the painting map, and the allele frequency model. Simulated genomes copy haplotypes only from the held-out people, with fixed seeds. The benchmark phases the same genome once with Beagle and sends those phased genotypes to both painters. It compares og projection with ADMIXTURE `-P` at 10,567 model SNVs, with bootstrap disabled.

## Regression check

Run the benchmark with the default chromosomes and five genomes per scenario.

```bash
THREADS=4 OG=<og-linux-binary> tests/bench/ancestry-engines.sh <panel_dir> <work_dir>
```

The gate allows a decrease of at most 0.01 from the committed og accuracy and dosage R² in each scenario. S0 also requires accuracy of at least 0.98. The table reports RFMix alongside og for comparison. Use a fresh work directory when you change the panel, binary, or chromosomes. Reusing a work directory keeps finished runs and scores their outputs again.

The measured run passed all five rows. With `BENCH_SHIFT=1`, scoring advances each og label by one group and every row fails. A second control rotated only og's two labels on chromosome 9 of one S3 genome. S3 failed the gate while the other four rows passed.

## Exploratory run with 400 references per group

Repainting the same phased inputs with og at a cap of 400 raised S2 accuracy from 0.8764 to 0.9060 and dosage R² from 0.8829 to 0.9081. S4 accuracy rose from 0.7817 to 0.8070 while dosage R² fell from 0.8049 to 0.7822. The exploratory run changed og's reference count alone. The equal-reference table and gate use the shipping cap of 200 for both engines.
