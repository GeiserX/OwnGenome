# Ancestry engine benchmark: og against RFMix and ADMIXTURE

> Result, 2026-11-24: **gate: FAIL.** og's painter loses to RFMix in all five scenarios. og's projector matches ADMIXTURE `-P` in every scenario. Step 26b should not move from RFMix to og's painter on this result.

Step 26b paints each stretch of a genome with a continental group (AFR, WEA, CSA, EAS, AMR) and gives the genome's continental fractions. Until now RFMix painted and ADMIXTURE `-P` gave the fractions. Our own engines, `og ancestry-fine paint` and `og ancestry-fine project`, are meant to replace both. A real genome has no known truth, so agreement with RFMix alone says nothing about which engine is right. This bench builds genomes whose truth is known and scores both engines against that truth.

## What the bench does

`tests/bench/ancestry-engines.sh <panel_dir> <work_dir>` runs everything. The comment at the top of the script describes every step and every column. In short:

1. It holds out 4 painting references per continental group, 20 people in all. It builds a bench model without them: the panel VCFs, each group's allele frequencies (`og ancestry-fine freqs --exclude`) and the painting map.
2. It simulates 5 genomes per scenario from the held-out people only, with fixed seeds (`og ancestry-fine simulate`). T is ancestry events per Morgan, about the generations since the mix.

   | Scenario | Groups | T |
   |---|---|---|
   | S0 | one group, 100% (one genome each for AFR, WEA, CSA, EAS, AMR) | none |
   | S1 | AFR 0.8, WEA 0.2 | 6 |
   | S2 | WEA 0.5, AMR 0.4, AFR 0.1 | 12 |
   | S3 | WEA 0.5, CSA 0.5 | 20 |
   | S4 | EAS 0.5, WEA 0.5 | 50 |

3. Beagle phases each genome against the bench panel. Then og and RFMix paint the same phased files, and `og ancestry-fine project` and ADMIXTURE `-P` give the fractions. The commands are step 26b's, with the bench model and no bootstrap.
4. It scores each painting against the truth, over every base of the truth's span, and checks the gate.

The genomes cover chr9, chr15 and chr21: 85,610 of the panel's SNVs and about 3.4 Morgans. The projection model keeps 10,567 SNVs on those chromosomes. Both painters use the same references: the panel's painting map without the held-out people, 46 per group (31 for AMR, which has no more).

**The gate**, per scenario: og's bp accuracy and dosage r2 are each at least RFMix's minus 0.01, and in S0 og's bp accuracy is at least 0.98. The script exits 1 when a row fails.

## Result

Five genomes per scenario, pooled. acc is the bp accuracy of the two copies (phase-free: both copies right, in either order). r2 is the dosage r2 of each group's copy count.

| Scenario | acc og | acc RFMix | r2 og | r2 RFMix | Gate |
|---|---|---|---|---|---|
| S0 | 0.6925 | 0.8985 | 0.7050 | 0.9028 | FAIL: acc, r2, S0 floor |
| S1 | 0.8137 | 0.9114 | 0.8486 | 0.9331 | FAIL: acc, r2 |
| S2 | 0.6743 | 0.8020 | 0.6498 | 0.8011 | FAIL: acc, r2 |
| S3 | 0.4023 | 0.5596 | 0.3570 | 0.5628 | FAIL: acc, r2 |
| S4 | 0.5529 | 0.6015 | 0.5102 | 0.6350 | FAIL: acc, r2 |

The script printed `GATE: FAIL` and exited 1.

The fractions pass. We measured the mean absolute error of the continental fractions against the truth. It is the same for `og ancestry-fine project` and ADMIXTURE `-P` in every scenario, to three decimals: 0.023, 0.017, 0.028, 0.046 and 0.034 for S0 to S4. The projector can replace ADMIXTURE `-P`.

## Why og's painter loses

| Scenario | hap og / RFMix | stretches truth / og / RFMix | missed og / RFMix / true | false og / RFMix | T simulated | T fitted |
|---|---|---|---|---|---|---|
| S0 | 0.8093 / 0.9391 | 6 / 259 / 16 | 0 / 0 / 0 | 1247 / 47 | none | 89.2 |
| S1 | 0.8304 / 0.8944 | 23 / 242 / 30 | 5 / 32 / 87 | 974 / 48 | 6 | 92.4 |
| S2 | 0.6478 / 0.7201 | 60 / 312 / 63 | 16 / 98 / 273 | 993 / 99 | 12 | 70.1 |
| S3 | 0.5448 / 0.6814 | 73 / 556 / 30 | 22 / 291 / 338 | 1849 / 70 | 20 | 132.6 |
| S4 | 0.6033 / 0.6857 | 189 / 552 / 57 | 43 / 619 / 918 | 829 / 48 | 50 | 132.3 |

Stretches are per genome. A breakpoint is where the pair of groups changes. A true breakpoint is missed when no call lies within 1 cM. A called breakpoint is false when no true one lies within 1 cM.

- og switches far too often. In the mixed scenarios it paints 3 to 11 times as many stretches as the truth has, and 43 times as many in S0. It calls 10 to 27 times as many false breakpoints as RFMix. RFMix errs the other way: once T reaches 20 it misses most true breakpoints (291 of 338 in S3, 619 of 918 in S4), while og misses few. Over the bases, og's extra switches cost more than RFMix's missed ones.
- The fitted T runs away. EM starts T at 10 and fits it between 70 and 133, whatever the simulated T is (none in S0, 50 in S4). The run-away is a symptom: in a first test on one S1 genome, holding T alone barely moved the accuracy. Holding T and the proportions together helps (next section).
- Phase is not the cause. hap is the per-haplotype accuracy after the better of the two phase alignments per chromosome. It is above og's phase-free accuracy in S0, S1, S3 and S4. Only in S2 is acc above hap by more than 0.02 (0.0265), and RFMix's own gap there is larger (0.0819). So a phase-flip layer would not close the loss. og loses on phase-free accuracy itself.
- RFMix does not reach S0's floor of 0.98 either. It reaches 0.8985. Per genome, og and RFMix score 0.9462 and 1.0000 (AFR), 0.6214 and 0.9809 (WEA), 0.3448 and 0.6667 (CSA), 0.7392 and 0.9855 (EAS), 0.8111 and 0.8594 (AMR). Some held-out people come from populations with known mixed ancestry: Bedouin for WEA, and Colombian, Maya and Pima for AMR. For them the "100% one group" truth is not exact. The floor stays as it is. og must reach it, and this bench shows how far both engines are from it.
- og is much cheaper. It paints a genome in 6 to 12 s with 67 to 139 MB at peak. RFMix takes 642 to 804 s with about 930 MB.

## What closes part of the gap

We repainted the same phased genomes with a diagnostic og build that is not part of this change. Every other input stayed the same, and RFMix's results did not change. Each cell is acc / r2, and a bold cell passes its condition.

| Variant | S0 | S1 | S2 | S3 | S4 |
|---|---|---|---|---|---|
| og as shipped | 0.6925 / 0.7050 | 0.8137 / 0.8486 | 0.6743 / 0.6498 | 0.4023 / 0.3570 | 0.5529 / 0.5102 |
| mu and T held at their start | 0.7702 / 0.7768 | 0.8886 / 0.9062 | 0.7143 / 0.6982 | 0.4814 / 0.4282 | **0.6278** / 0.6199 |
| 400 painting references per group | 0.8020 / 0.8392 | 0.8984 / 0.9221 | 0.7452 / 0.7371 | **0.5966 / 0.5980** | **0.7377 / 0.7585** |
| both | 0.8575 / 0.8800 | **0.9275 / 0.9446** | 0.7548 / 0.7425 | **0.6682 / 0.6635** | **0.7729 / 0.7901** |
| RFMix | 0.8985 / 0.9028 | 0.9114 / 0.9331 | 0.8020 / 0.8011 | 0.5596 / 0.5628 | 0.6015 / 0.6350 |

"mu and T held" keeps the proportions at the continental projection and T at its start of 10 through every EM round. EM does not fit them again from the posterior. "400 painting references" gives og up to 400 references per group (31 for AMR), and RFMix keeps its 46. That comparison is not like for like, because RFMix might also gain from more references. It is still the comparison that matters for step 26b. With 400 references per group og paints a genome in about 45 s, about 15 times faster than RFMix with 46.

Both changes together pass S1, S3 and S4 and still fail S0 and S2. S0 reaches 0.8575 against the 0.98 floor, and og still paints 46 stretches per genome where the truth has 6. We have not yet found the cause of the rest of the gap. Two more tries in the diagnostic build did not help S0:

- Each ancestry copied haplotypes of the other groups 5% of the time, a fixed copy matrix as FLARE uses. The S0 AFR genome went from 0.9462 to 0.9284.
- A model on the groups' allele frequencies alone scored 0.3319 in S0.

## Controls

A gate must be able to fail, and it must be able to pass.

- Red control, `BENCH_SHIFT=1`: og's labels move one group along before scoring. The run exits 1 and every row fails, for example `S0 0.0110 0.8985 0.0314 0.9028 ... FAIL: acc r2 S0-acc<0.98`.
- Green control, `BENCH_SELF=1`: RFMix's painting is scored in og's place. Both columns are equal, so the acc and r2 conditions pass in every scenario. The run exits 1 from the S0 floor only: `S0 0.8985 0.8985 0.9028 0.9028 ... FAIL: S0-acc<0.98`.

## How it ran

- Code: OwnGenome at `05bd7cfb`, the head of the painter wiring change at the time. The engine (`rust/crates/og-admix`) and the painting maths are the same at `a86899e1`. og is built with `cargo zigbuild --release --locked -p og --target x86_64-unknown-linux-musl`.
- Panel: gnomAD's HGDP + 1000 Genomes panel (CC0) as `scripts/setup.sh --ancestry-fine-panel` installs it. The bench reads no personal genome, so the genomes and the numbers can be published.
- Images, from `versions.env`: `quay.io/biocontainers/rfmix:2.03.r0.9505bfa--h503566f_8`, `quay.io/biocontainers/beagle:5.5_27Feb25.75f--hdfd78af_0`, `quay.io/biocontainers/admixture:1.3.0--0`, `pgscatalog/plink2:2.00a5.10`, `staphb/bcftools:1.24`, and og in `ubuntu:24.04`.
- Machine: Intel Core i5-14500, each tool limited to 4 threads (`THREADS=4`, `docker run --cpus 4`). The 25 genomes took about 5.5 hours, nearly all of it in RFMix.
- Run: `THREADS=4 OG=<og binary> tests/bench/ancestry-engines.sh <panel_dir> <work_dir>`. A second run keeps every finished step and only scores again, so each control takes a few minutes.
