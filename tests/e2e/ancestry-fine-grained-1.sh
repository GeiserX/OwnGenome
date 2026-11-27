#!/usr/bin/env bash
# Step 26b (fine-grained ancestry) against a reduced panel, read from gnomAD's
# bucket by setup.sh --ancestry-fine-panel for the fixture's slices plus
# chr22:20-35 Mb (the slices are too short to paint: the painting needs 5 cM):
#   - without the panel, step 26b says so in one line and exits 0;
#   - the panel installs: both levels of the model, painting references, and
#     none of the three people left out (A, B and C below) in it. Oceania, which has one
#     unrelated person in the panel's metadata, is left out of the model at
#     both levels (5 continental and 17 regional groups), and the install log
#     and HG002T's report name it;
#   - HG002 from its gVCF on the slices (cases 21 and 03): the continental
#     fractions sum to 1, each carries a bootstrap range, and with no
#     chromosome long enough the step says it painted nothing. Only about a
#     quarter of the panel's SNVs are in the slices, too few to name a group;
#   - HG002T, HG002's GIAB v4.2.1 genotypes at the panel's SNVs inside GIAB's
#     confident regions (0/0 where GIAB lists no variant): West Eurasia first
#     with 70% or more, and the largest share of the painted chr22, half or more;
#     og painted chr22 with no failure and wrote its fitted parameters;
#   - MIX, a planted person: chr22 below 27.5 Mb and the first half of each
#     slice from a Yoruba (NA18486), the rest from a Han Chinese (NA18525),
#     both left out of the panel: Africa and East Asia each 25% to 75%, the
#     painting African left of 25.5 Mb and East Asian right of 29.5 Mb, and
#     neither side passing the other group's check;
#   - SIM, three genomes that og ancestry-fine simulate builds from A, B and C
#     (NA06984, a Utah resident of northern and western European ancestry,
#     also left out of the panel), with their truth: Africa and East Asia at T 6, West
#     Eurasia and East Asia at T 20, and one West Eurasian with no mix. On
#     each, the painting of chr22:20-35 Mb must have the copies of each group
#     right on enough of its bases and few false breakpoints per Morgan, and
#     the fractions must be near the truth (floors measured on this panel,
#     with a margin); a painting with its breakpoints moved must fail them;
#   - TWO, a two-sample gVCF whose first sample is MIX and second HG002T: the
#     step reads the first sample column, so it must stop with the reason and
#     write no fractions, instead of reporting MIX's ancestry under TWO's name;
#   - red control: the same panel with every reference person's labels moved
#     to the next group (Africa's people called West Eurasia, and so on) and
#     the group frequencies counted again: on it the HG002T fraction checks
#     fail (sample HG002R, the same genotypes; the red copy then loses the
#     genetic map of chr22, so that run only estimates fractions), and before
#     that SIMAE against it (sample SIMAR) fails SIMAE's floors;
#   - the Nextflow module on HG002T (with html_report) writes the fractions and
#     the painting of the bash step, and the HTML report has its card; on
#     TWO's gVCF it stops the same way and publishes nothing.
# The time and disk of the panel install, and the scores and time of the
# simulated genomes, go to the job summary.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
G="$GENOME_DIR"
PANEL="${G}/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}"
A=NA18486 B=NA18525 C=NA06984
LONG_CHR=chr22 LONG_START=20000000 LONG_END=35000000 JUNCTION=27500000

frac() {  # frac SAMPLE LEVEL CODE: the fraction column of that row
  awk -F'\t' -v l="$2" -v c="$3" '$1 == l && $2 == c {print $4}' "${G}/$1/ancestry_fine/$1_ancestry_fine.tsv" 2>/dev/null
}
ge() { awk -v x="${1:-}" -v y="$2" 'BEGIN {exit !(x ~ /^[0-9.]+$/ && x + 0 >= y)}'; }
between() { awk -v x="${1:-}" -v lo="$2" -v hi="$3" 'BEGIN {exit !(x ~ /^[0-9.]+$/ && x + 0 >= lo && x + 0 <= hi)}'; }
# side SAMPLE CHROM FROM TO CODE: share of the painted bp of CHROM between FROM and TO painted CODE
side() {
  awk -F'\t' -v c="$2" -v a="$3" -v b="$4" -v g="$5" '
    $1 == c { s = ($2 > a ? $2 : a); e = ($3 < b ? $3 : b); if (e > s) { t += e - s; if ($5 == g) h += e - s } }
    END { if (t > 0) printf "%.3f\n", h / t }' "${G}/$1/ancestry_fine/$1_painting.bed" 2>/dev/null
}

# --- no panel: one line, exit 0 ---------------------------------------------------
ANCESTRY_FINE_DIR=none run_step 26b-ancestry-fine.sh "$SAMPLE"
check_step_exit "26b-ancestry-fine.sh (no panel)"
check_eq "step 26b without a panel prints one line" "$(grep -c . "$STEP_LOG")" 1
check "that line says it was skipped and how to install the panel" has '^Step 26b skipped: .*setup\.sh --ancestry-fine-panel' "$(cat "$STEP_LOG")"

# --- the reduced panel ----------------------------------------------------------------
REGIONS="${CASE_TMP}/regions.bed"
{ cat "${FIXTURE_DIR}/regions.bed"; printf '%s\t%s\t%s\n' "$LONG_CHR" "$LONG_START" "$LONG_END"; } > "$REGIONS"
echo "+ ANCESTRY_FINE_REGIONS=... ANCESTRY_FINE_EXCLUDE=${A},${B},${C} scripts/setup.sh --ancestry-fine-panel"
T0=$(date +%s)
ANCESTRY_FINE_REGIONS="$REGIONS" ANCESTRY_FINE_EXCLUDE="${A},${B},${C}" \
  "${REPO}/scripts/setup.sh" --ancestry-fine-panel "$G" 2>&1 | tee "$STEP_LOG"
check_eq "setup.sh --ancestry-fine-panel exits 0" "${PIPESTATUS[0]}" 0
T_SETUP=$(( $(date +%s) - T0 ))
check "the panel is finished (model/DONE)" test -s "${PANEL}/model/DONE"
check_eq "continental level: 5 groups (Oceania left out)" "$(awk -F'\t' '$1 == "continent" {print $2}' "${PANEL}/model/levels.tsv" 2>/dev/null)" 5
check_eq "regional level: 17 groups (Oceania left out)" "$(awk -F'\t' '$1 == "region" {print $2}' "${PANEL}/model/levels.tsv" 2>/dev/null)" 17
check "Oceania has no reference person and no painting reference" \
  awk -F'\t' '$3 == "OCE" && ($5 == 1 || $6 == 1) {bad = 1} END {exit bad}' "${PANEL}/samples.tsv"
check "the install log names Oceania as left out" \
  has 'Left out of the model, fewer than 10 unrelated reference people: Oceania \(continental and regional group\)' "$(cat "$STEP_LOG")"
check_ge "SNVs of the model" "$(wc -l < "${PANEL}/model/ref.bim" 2>/dev/null | tr -d ' ')" 500
check_ge "painting references" "$(wc -l < "${PANEL}/model/paint_map.tsv" 2>/dev/null | tr -d ' ')" 150
# shellcheck disable=SC2016 # $1 is awk's field
check "the three people left out (the mix and the simulated genomes) are not in the panel" \
  awk -F'\t' -v a="$A" -v b="$B" -v c="$C" '$1 == a || $1 == b || $1 == c {bad = 1} END {exit bad}' "${PANEL}/samples.tsv"
printf '#### Step 26b reduced panel (fixture slices plus %s:%s-%s)\n\n%s s to install, %s on disk, %s SNVs, %s in the model\n\n' \
  "$LONG_CHR" "$LONG_START" "$LONG_END" "$T_SETUP" "$(du -sh "$PANEL" | cut -f1)" \
  "$(wc -l < "${PANEL}/panel/sites.tsv" | tr -d ' ')" "$(wc -l < "${PANEL}/model/ref.bim" | tr -d ' ')" >> "$E2E_NOTES"

# --- HG002 from its gVCF on the fixture's slices ---------------------------------------
T0=$(date +%s)
run_step 26b-ancestry-fine.sh "$SAMPLE"
check_step_exit "26b-ancestry-fine.sh ${SAMPLE}"
echo "step 26b on ${SAMPLE}: $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
cat "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.tsv" 2>/dev/null
check "it read the gVCF" has 'Input: the gvcf' "$(cat "$STEP_LOG")"
SUM=$(awk -F'\t' '$1 == "continent" {s += $4} END {printf "%.3f", s}' "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.tsv" 2>/dev/null)
check "the continental fractions sum to 1 (${SUM})" between "$SUM" 0.99 1.01
# Only about a quarter of the model's SNVs fall in the fixture's
# slices and are called in the gVCF, too few to name a group: this run checks
# the shape of the output, and HG002T below checks where the genome lands.
check_ge "SNVs of the panel genotyped from the gVCF" \
  "$(awk '{print $1; exit}' <<< "$(grep -o '[0-9]* of the panel' "$STEP_LOG")")" 500
check "each continental fraction has a bootstrap range" \
  awk -F'\t' '$1 == "continent" && ($6 == "." || $7 == "." || $6 > $4 + 1e-9 || $7 < $4 - 1e-9) {bad = 1} END {exit bad}' \
  "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.tsv"
check "no slice spans 5 cM, so nothing was painted, and the summary says so" \
  has '^Chromosome painting: none' "$(cat "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_ancestry_fine.txt" 2>/dev/null)"
check "the PNG is a PNG" has '^89504e47' "$(head -c 4 "${G}/${SAMPLE}/ancestry_fine/${SAMPLE}_painting.png" 2>/dev/null | od -An -tx1 | tr -d ' \n')"

# --- HG002T: GIAB's genotypes at the panel's SNVs ------------------------------------
# Inside GIAB's confident regions a site GIAB lists no variant at is 0/0.
mkdir -p "${G}/HG002T/vcf"
docker run --rm -u "$(id -u):$(id -g)" -v "${G}:/genome:ro" -v "${FIXTURE_DIR}:/fx:ro" "$BCFTOOLS_IMAGE" \
  bcftools view -R "/genome/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/panel/positions.tsv" -Ov \
  /fx/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > "${CASE_TMP}/truth_at_sites.vcf"
"${REPO}/bin/og" dev fixture truth-genotypes "${PANEL}/panel/sites.tsv" "${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai" \
  "HG002T=${CASE_TMP}/truth_at_sites.vcf,${FIXTURE_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
  > "${G}/HG002T/vcf/HG002T.vcf"
in_genome "$BCFTOOLS_IMAGE" bash -c 'bcftools view -Oz -o HG002T/vcf/HG002T.vcf.gz HG002T/vcf/HG002T.vcf && bcftools index -t HG002T/vcf/HG002T.vcf.gz'
check_ge "HG002T: GIAB genotypes at panel SNVs" "$(vcf_count HG002T/vcf/HG002T.vcf.gz)" 3000
T0=$(date +%s)
run_step 26b-ancestry-fine.sh HG002T
check_step_exit "26b-ancestry-fine.sh HG002T"
echo "step 26b on HG002T (painting ${LONG_CHR}): $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
cat "${G}/HG002T/ancestry_fine/HG002T_ancestry_fine.txt" 2>/dev/null
WEA_T=$(frac HG002T continent WEA) PAINT_T=$(side HG002T "$LONG_CHR" "$LONG_START" "$LONG_END" WEA)
check_eq "HG002T: West Eurasia comes first" "$(awk -F'\t' '$1 == "continent" {print $2; exit}' "${G}/HG002T/ancestry_fine/HG002T_ancestry_fine.tsv" 2>/dev/null)" WEA
check "HG002T: with 70% or more (${WEA_T})" ge "$WEA_T" 0.7
check_eq "HG002T: West Eurasia is the largest share of the painting" \
  "$(awk -F'\t' '$1 == "painting" {print $2; exit}' "${G}/HG002T/ancestry_fine/HG002T_ancestry_fine.tsv" 2>/dev/null)" WEA
check "HG002T: half or more of ${LONG_CHR} painted West Eurasia (${PAINT_T:-none})" ge "$PAINT_T" 0.5
check_ge "HG002T: painted Mb of ${LONG_CHR} on two haplotypes" \
  "$(awk -F'\t' -v c="$LONG_CHR" '$1 == c {t += $3 - $2} END {print int(t / 1e6)}' "${G}/HG002T/ancestry_fine/HG002T_painting.bed" 2>/dev/null)" 20
check "HG002T: the report says Oceania was left out of the model" \
  has '^Left out of the model: Oceania \(continental and regional group\)' "$(cat "${G}/HG002T/ancestry_fine/HG002T_ancestry_fine.txt" 2>/dev/null)"
check "HG002T: and the table has its left_out rows" \
  awk -F'\t' '$1 == "left_out" && $2 == "OCE" {n++} END {exit n != 2}' "${G}/HG002T/ancestry_fine/HG002T_ancestry_fine.tsv"
# og ancestry-fine project wrote each level's fit, after a bootstrap that drew
# the reference people again.
for l in continent region; do
  # shellcheck disable=SC2016 # $1 and $2 are awk's fields
  check "HG002T: the ${l} fit drew 100 replicates with the reference people drawn again" \
    awk -F'\t' '$1 == "replicates" {r = $2} $1 == "reference_resampling" {s = $2} END {exit !(r == 100 && s == "yes")}' \
    "${G}/HG002T/ancestry_fine/work/fractions/HG002T_${l}.fit.tsv"
done
# og ancestry-fine paint painted every phased chromosome with our own model.
check "HG002T: og painted ${LONG_CHR} with no failure" test ! -e "${G}/HG002T/ancestry_fine/work/paint/failed.txt"
# shellcheck disable=SC2016 # $1 and $2 are awk's fields
check "HG002T: the painting wrote its fitted parameters" \
  awk -F'\t' '$1 == "em_rounds" {r = 1} $1 == "mu" && $2 == "WEA" {w = 1} END {exit !(r && w)}' \
  "${G}/HG002T/ancestry_fine/work/paint/params.tsv"

# --- MIX: half a Yoruba, half a Han Chinese ---------------------------------------------
# Their phased genotypes at the panel's SNVs, and those of C for the
# simulated genomes below, from the same bucket (the panel holds none of
# them), kept in MIX/tmp/held/<chrom>.vcf.gz. MIX gets A's and B's unphased:
# the step phases them itself.
mkdir -p "${G}/MIX/vcf" "${G}/MIX/tmp/held"
# shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
docker run --rm -u "$(id -u):$(id -g)" -v "${G}:/genome" -w /genome/MIX/tmp "$BCFTOOLS_HTTPS_IMAGE" bash -c '
  for c in $(cut -f1 "$1" | uniq); do
    r=$(awk -v c="$c" '"'"'$1 == c {printf "%s%s:%d-%d", s, $1, $2 + 1, $3; s = ","}'"'"' "$3")
    bcftools view -s "$2" -T "$1" -r "$r" -m2 -M2 -v snps \
      "https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_${c}.filtered.SNV_INDEL.phased.shapeit5.bcf" \
      -Oz -o "held/${c}.vcf.gz" 2>/dev/null && bcftools index -t "held/${c}.vcf.gz" \
      && bcftools query -H -s "$4" -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n" "held/${c}.vcf.gz"
  done; rm -f ./*.csi' _ "/genome/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/panel/positions.tsv" "${A},${B},${C}" \
  "/genome/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/regions.bed" "${A},${B}" > "${CASE_TMP}/ab.tsv"
check_ge "genotypes of ${A} and ${B} at the panel's SNVs" "$(grep -vc '^#' "${CASE_TMP}/ab.tsv" || true)" 3000
# The first half of each region (of chr22:20-35 Mb, below 27.5 Mb) from A,
# the rest from B; the columns are found by name in bcftools' -H header.
{ printf '##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n'
  cut -f1 "${PANEL}/panel/sites.tsv" | uniq | awk -F'\t' 'NR == FNR {u[$1] = 1; next} ($1 in u) {printf "##contig=<ID=%s,length=%s>\n", $1, $2}' - "${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai"
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMIX\n'
  awk -F'\t' -v OFS='\t' -v regions="$REGIONS" -v a="$A" -v b="$B" '
    BEGIN { while ((getline l < regions) > 0) { split(l, f, "\t"); n++; rc[n] = f[1]; rs[n] = f[2]; re[n] = f[3] } }
    /^#/ { for (i = 5; i <= NF; i++) { if (index($i, "]" a ":")) ca = i; if (index($i, "]" b ":")) cb = i }; next }
    {
      mid = 0
      for (i = 1; i <= n; i++) if (rc[i] == $1 && $2 > rs[i] && $2 <= re[i]) mid = (rs[i] + re[i]) / 2
      g = ($2 < mid) ? $ca : $cb; sub(/\|/, "/", g); if (g == "1/0") g = "0/1"
      print $1, $2, ".", $3, $4, ".", "PASS", ".", "GT", g
    }' "${CASE_TMP}/ab.tsv"; } > "${G}/MIX/vcf/MIX.vcf"
in_genome "$BCFTOOLS_IMAGE" bash -c 'bcftools sort -Oz -o MIX/vcf/MIX.vcf.gz MIX/vcf/MIX.vcf && bcftools index -t MIX/vcf/MIX.vcf.gz'
T0=$(date +%s)
run_step 26b-ancestry-fine.sh MIX
check_step_exit "26b-ancestry-fine.sh MIX"
echo "step 26b on MIX: $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
cat "${G}/MIX/ancestry_fine/MIX_ancestry_fine.txt" 2>/dev/null
check "MIX: Africa 25% to 75% ($(frac MIX continent AFR))" between "$(frac MIX continent AFR)" 0.25 0.75
check "MIX: East Asia 25% to 75% ($(frac MIX continent EAS))" between "$(frac MIX continent EAS)" 0.25 0.75
L=$(side MIX "$LONG_CHR" "$LONG_START" $((JUNCTION - 2000000)) AFR)
R=$(side MIX "$LONG_CHR" $((JUNCTION + 2000000)) "$LONG_END" EAS)
check "MIX: Africa on 70% or more of ${LONG_CHR} left of the junction (${L:-none})" ge "$L" 0.7
check "MIX: East Asia on 60% or more of ${LONG_CHR} right of the junction (${R:-none})" ge "$R" 0.6
# The two sides swapped: the same checks on the other group must not pass, so
# neither of them is a check that cannot fail.
check "MIX: the left side is not East Asian ($(side MIX "$LONG_CHR" "$LONG_START" $((JUNCTION - 2000000)) EAS))" \
  test "$(ge "$(side MIX "$LONG_CHR" "$LONG_START" $((JUNCTION - 2000000)) EAS)" 0.7 && echo pass)" != pass
check "MIX: the right side is not African ($(side MIX "$LONG_CHR" $((JUNCTION + 2000000)) "$LONG_END" AFR))" \
  test "$(ge "$(side MIX "$LONG_CHR" $((JUNCTION + 2000000)) "$LONG_END" AFR)" 0.6 && echo pass)" != pass

# --- SIM: simulated genomes with known truth ------------------------------------------
# og ancestry-fine simulate builds each genome from the phased haplotypes of
# A, B and C (fetched above, none of them in the panel) along the genetic map,
# and writes the truth of each haplotype beside it. Fixed seeds, so the
# genomes, their phasing and their painting are the same on every run. The
# step runs on each as on any person (Beagle phases it first). On LONG_CHR,
# the one chromosome long enough to paint, the painting is scored against
# the truth: the share of the bases whose copies of each group are right,
# and the false breakpoints per Morgan; the continental fractions are scored
# against the truth's share of the bases of every simulated chromosome. The
# floors and ceilings come from the scores measured on this reduced panel,
# with a margin: accuracy 0.03 below the measured share and fraction error
# 0.03 above; false breakpoints one more than measured, per Morgan of the
# painted region (0.29 Morgan, so one breakpoint is 3.4 per Morgan), rounded
# up. Measured: SIMAE 0.9837, 0 false of 1 called, 0.0867; SIMWE 0.9557, 2
# of 5, 0.0466; SIMW 1.0000, 0 of 0, 0.0645. Two red controls must fail
# them: SIMAE against the rotated labels (in the red control below), and
# SIMWE's painting with every breakpoint moved 1.5 Mb (about 3 cM).
#   name   mu                T   seed  acc   false/M  |q-t|
SIM_CASES="SIMAE  AFR=0.5,EAS=0.5  6   1     0.9537  3.5      0.1167
SIMWE  WEA=0.5,EAS=0.5  20  2     0.9257  10.3     0.0766
SIMW   WEA=1            1   3     0.9700  3.5      0.0945"
printf '%s\tAFR\n%s\tEAS\n%s\tWEA\n' "$A" "$B" "$C" > "${CASE_TMP}/held_people.tsv"
# sim_score SAMPLE TRUTH [SHIFT]: SAMPLE's painting of LONG_CHR from
# LONG_START to LONG_END against TRUTH, as "accuracy false-per-Morgan false
# called true Morgans". A base is right when its two painted groups are its
# two true groups in either order; a base painted on neither haplotype is
# wrong. A breakpoint is a base where that pair of groups changes; a called
# one with no true one within 1 cM is false. SHIFT moves the whole painting
# that many bases to the right (a red control).
sim_score() {
  awk -F'\t' -v c="$LONG_CHR" -v lo="$LONG_START" -v hi="$LONG_END" -v sh="${3:-0}" -v map="${PANEL}/maps/${LONG_CHR}.map" '
    function cm(p,   l, h, m) {
      if (p <= mp[1]) return mc[1]
      if (p >= mp[nm]) return mc[nm]
      l = 1; h = nm
      while (h - l > 1) { m = int((l + h) / 2); if (mp[m] <= p) l = m; else h = m }
      return mp[h] == mp[l] ? mc[l] : mc[l] + (mc[h] - mc[l]) * (p - mp[l]) / (mp[h] - mp[l])
    }
    function grp(k, h, x,   i) {
      for (i = 1; i <= n[k, h]; i++) if (s[k, h, i] <= x && x < e[k, h, i]) return g[k, h, i]
      return "-"
    }
    function pair(x, y) { return x < y ? x "/" y : y "/" x }
    BEGIN { while ((getline l < map) > 0) { split(l, f, /[ \t]+/); if (f[4] != "") { nm++; mp[nm] = f[4] + 0; mc[nm] = f[3] + 0 } } }
    FNR == 1 { file++ }
    /^#/ || $1 != c { next }
    {
      k = file == 1 ? "t" : "c"; h = file == 1 ? $4 + 1 : $4
      a = $2 + (k == "c" ? sh : 0); b = $3 + (k == "c" ? sh : 0)
      if (a < lo) a = lo
      if (b > hi) b = hi
      if (b <= a) next
      i = ++n[k, h]; s[k, h, i] = a; e[k, h, i] = b; g[k, h, i] = $5
      if (!(a in seen)) { seen[a] = 1; p[++np] = a }
      if (!(b in seen)) { seen[b] = 1; p[++np] = b }
    }
    END {
      for (i = 2; i <= np; i++) { v = p[i]; for (j = i - 1; j >= 1 && p[j] > v; j--) p[j + 1] = p[j]; p[j + 1] = v }
      for (i = 1; i < np; i++) {
        x = p[i]; t1 = grp("t", 1, x); t2 = grp("t", 2, x)
        if (t1 == "-" || t2 == "-") continue
        c1 = grp("c", 1, x); c2 = grp("c", 2, x)
        tp = pair(t1, t2); cp = (c1 == "-" || c2 == "-") ? "" : pair(c1, c2)
        if (first == "") first = x
        last = p[i + 1]; tot += p[i + 1] - x
        if (cp == tp) ok += p[i + 1] - x
        if (ltp != "" && tp != ltp) tb[++nt] = cm(x)
        if (lcp != "" && cp != "" && cp != lcp) cb[++nc] = cm(x)
        ltp = tp; if (cp != "") lcp = cp
      }
      for (i = 1; i <= nc; i++) { hit = 0; for (j = 1; j <= nt; j++) if ((cb[i] - tb[j]) ^ 2 <= 1) hit = 1; if (!hit) nf++ }
      mo = (tot > 0) ? (cm(last) - cm(first)) / 100 : 0
      acc = (tot > 0) ? ok / tot : 0; fpm = (mo > 0) ? nf / mo : 0
      printf "%.4f %.2f %d %d %d %.3f\n", acc, fpm, nf + 0, nc + 0, nt + 0, mo
    }' "$2" "${G}/$1/ancestry_fine/$1_painting.bed" 2>/dev/null
}
# sim_qerr SAMPLE TRUTH: the largest |continental fraction - the truth's share
# of the bases| over the groups of either.
sim_qerr() {
  awk -F'\t' 'FNR == 1 { file++ }
    file == 1 && !/^#/ { bp[$5] += $3 - $2; tot += $3 - $2; next }
    file == 2 && $1 == "continent" { q[$2] = $4 }
    END {
      if (tot == 0) exit
      for (x in q) { d = q[x] - (x in bp ? bp[x] / tot : 0); if (d < 0) d = -d; if (d > m) m = d }
      for (x in bp) if (!(x in q)) m = 1
      printf "%.4f\n", m
    }' "$2" "${G}/$1/ancestry_fine/$1_ancestry_fine.tsv" 2>/dev/null
}
# below / above VALUE LIMIT: VALUE is a number under / over LIMIT.
below() { awk -v x="${1:-}" -v y="$2" 'BEGIN {exit !(x ~ /^[0-9.]+$/ && x + 0 < y)}'; }
T_SIM=$(date +%s)
printf '#### Step 26b on simulated genomes (%s:%s-%s painted)\n\n| genome | accuracy (floor) | false breakpoints per Morgan (ceiling) | max abs fraction error (ceiling) | breakpoints called / true |\n|---|---|---|---|---|\n' \
  "$LONG_CHR" "$LONG_START" "$LONG_END" >> "$E2E_NOTES"
# The loop reads its rows on fd 3: docker run -i in in_genome reads stdin.
while read -r n mu t seed acc_min fbpm_max q_max <&3; do
  mkdir -p "${G}/${n}/vcf"
  echo "+ og ancestry-fine simulate --mu ${mu} --generations ${t} --seed ${seed} --sample ${n}"
  "${REPO}/bin/og" ancestry-fine simulate --panel "${G}/MIX/tmp/held" --maps "${PANEL}/maps" \
    --people "${CASE_TMP}/held_people.tsv" --mu "$mu" --generations "$t" --seed "$seed" --sample "$n" \
    --out "${G}/${n}/vcf/${n}.vcf" --truth "${CASE_TMP}/${n}.truth.bed"
  check_eq "${n}: og ancestry-fine simulate exits 0" "$?" 0
  # shellcheck disable=SC2016  # $1 belongs to the inner bash
  in_genome "$BCFTOOLS_IMAGE" bash -c 'bcftools view -Oz -o "$1/vcf/$1.vcf.gz" "$1/vcf/$1.vcf" && bcftools index -t "$1/vcf/$1.vcf.gz"' _ "$n"
  ANCESTRY_FINE_BOOTSTRAP=0 run_step 26b-ancestry-fine.sh "$n"
  check_step_exit "26b-ancestry-fine.sh ${n}"
  check "${n}: og painted ${LONG_CHR} with no failure" test ! -e "${G}/${n}/ancestry_fine/work/paint/failed.txt"
  read -r s_acc s_fbpm _ s_called s_true _ <<< "$(sim_score "$n" "${CASE_TMP}/${n}.truth.bed")"
  s_q=$(sim_qerr "$n" "${CASE_TMP}/${n}.truth.bed")
  echo "${n}: accuracy ${s_acc:-none}, false breakpoints per Morgan ${s_fbpm:-none} (${s_called:-?} called, ${s_true:-?} true), max |q - truth| ${s_q:-none}"
  check "${n}: the copies of each group are right on ${acc_min} or more of the bases (${s_acc:-none})" ge "$s_acc" "$acc_min"
  check "${n}: at most ${fbpm_max} false breakpoints per Morgan (${s_fbpm:-none})" \
    test "$(ge "$s_fbpm" 0 && ! below "$fbpm_max" "$s_fbpm" && echo pass)" = pass
  check "${n}: each continental fraction within ${q_max} of the truth (${s_q:-none})" \
    test "$(ge "$s_q" 0 && ! below "$q_max" "$s_q" && echo pass)" = pass
  printf '| %s (%s, T %s) | %s (%s) | %s (%s) | %s (%s) | %s / %s |\n' "$n" "$mu" "$t" "$s_acc" "$acc_min" \
    "$s_fbpm" "$fbpm_max" "$s_q" "$q_max" "$s_called" "$s_true" >> "$E2E_NOTES"
done 3<<< "$SIM_CASES"
check_ge "SIMWE: the truth has breakpoints on ${LONG_CHR} to find" \
  "$(sim_score SIMWE "${CASE_TMP}/SIMWE.truth.bed" | awk '{print $5}')" 2
# Red control: SIMWE's painting moved 1.5 Mb to the right, so each called
# breakpoint is about 3 cM from where it belongs, fails its floors.
read -r r_acc r_fbpm _ <<< "$(sim_score SIMWE "${CASE_TMP}/SIMWE.truth.bed" 1500000)"
WE_ACC=$(awk '$1 == "SIMWE" {print $5}' <<< "$SIM_CASES") WE_FBPM=$(awk '$1 == "SIMWE" {print $6}' <<< "$SIM_CASES")
check "red: SIMWE's painting moved 1.5 Mb fails the accuracy floor (${r_acc:-none} < ${WE_ACC})" below "$r_acc" "$WE_ACC"
check "red: and the false-breakpoint ceiling (${r_fbpm:-none} > ${WE_FBPM})" \
  test "$(ge "$r_fbpm" 0 && below "$WE_FBPM" "$r_fbpm" && echo pass)" = pass
printf '\nSimulated genomes: %s s for the three\n\n' "$(( $(date +%s) - T_SIM ))" >> "$E2E_NOTES"

# --- TWO: a two-sample gVCF ------------------------------------------------------------
# Its VCF is HG002T's (one sample), so the check that stops the run is the gVCF's.
mkdir -p "${G}/TWO/vcf"
in_genome "$BCFTOOLS_IMAGE" bash -c 'bcftools merge -Oz -o TWO/vcf/TWO.g.vcf.gz MIX/vcf/MIX.vcf.gz HG002T/vcf/HG002T.vcf.gz && bcftools index -t TWO/vcf/TWO.g.vcf.gz'
cp "${G}/HG002T/vcf/HG002T.vcf.gz" "${G}/TWO/vcf/TWO.vcf.gz"
cp "${G}/HG002T/vcf/HG002T.vcf.gz.tbi" "${G}/TWO/vcf/TWO.vcf.gz.tbi"
check_eq "TWO: the gVCF holds MIX then HG002T" "$(bcf query -l TWO/vcf/TWO.g.vcf.gz | paste -sd, -)" "MIX,HG002T"
ANCESTRY_FINE_BOOTSTRAP=0 run_step 26b-ancestry-fine.sh TWO
check "TWO: step 26b stops on the two-sample gVCF (exit ${STEP_RC})" test "$STEP_RC" -ne 0
check "TWO: and says why" has 'holds 2 samples; step 26b reads a one-sample file' "$(cat "$STEP_LOG")"
check "TWO: and writes no fractions" test ! -e "${G}/TWO/ancestry_fine/TWO_ancestry_fine.tsv"

# --- red control: the reference labels moved to the next group ------------------------
RED="${G}/reference/ancestry_fine/red-control"
rm -rf "$RED"
cp -a "$PANEL" "$RED"
rm -rf "${RED}/model"
awk -F'\t' -v OFS='\t' '
  BEGIN { split("AFR WEA CSA EAS AMR OCE", c, " "); for (i = 1; i <= 6; i++) nc[c[i]] = c[i % 6 + 1]
          split("WAF EAF CAF SAF NWE NEE IBE ITA SAR CAU LEV BDN NAF CSA EAS SIB AMR OCE", r, " "); for (i = 1; i <= 18; i++) nr[r[i]] = r[i % 18 + 1] }
  NR == 1 { print; next }
  { if ($3 in nc) $3 = nc[$3]; if ($4 in nr) $4 = nr[$4]; print }' "${PANEL}/samples.tsv" > "${CASE_TMP}/rotated.tsv"
rm -f "${RED}/samples.tsv"
ANCESTRY_FINE_DIR="$RED" ANCESTRY_FINE_REGIONS="$REGIONS" ANCESTRY_FINE_SAMPLE_MAP="${CASE_TMP}/rotated.tsv" \
  "${REPO}/scripts/setup.sh" --ancestry-fine-panel "$G" > "${CASE_TMP}/red_setup.log" 2>&1
RED_RC=$?
# SIMAE against the rotated labels, while the red copy still has the map of
# the long chromosome: its painting and fractions must fail SIMAE's floors.
mkdir -p "${G}/SIMAR/vcf"
cp "${G}/SIMAE/vcf/SIMAE.vcf.gz" "${G}/SIMAR/vcf/SIMAR.vcf.gz"
cp "${G}/SIMAE/vcf/SIMAE.vcf.gz.tbi" "${G}/SIMAR/vcf/SIMAR.vcf.gz.tbi"
ANCESTRY_FINE_DIR="$RED" ANCESTRY_FINE_BOOTSTRAP=0 run_step 26b-ancestry-fine.sh SIMAR
check_step_exit "26b-ancestry-fine.sh SIMAR (SIMAE against the rotated labels)"
read -r r_acc _ <<< "$(sim_score SIMAR "${CASE_TMP}/SIMAE.truth.bed")"
r_q=$(sim_qerr SIMAR "${CASE_TMP}/SIMAE.truth.bed")
AE_ACC=$(awk '$1 == "SIMAE" {print $5}' <<< "$SIM_CASES") AE_Q=$(awk '$1 == "SIMAE" {print $7}' <<< "$SIM_CASES")
check "red: SIMAE against the rotated labels fails the accuracy floor (${r_acc:-none} < ${AE_ACC})" below "$r_acc" "$AE_ACC"
check "red: and the fraction ceiling (${r_q:-none} > ${AE_Q})" test "$(ge "$r_q" 0 && below "$AE_Q" "$r_q" && echo pass)" = pass
# The genetic map of the long chromosome is then taken out of the red copy, so
# no chromosome passes the 5 cM rule and the red run stops after the fractions:
# this control is for the fraction checks, and the painting checks have theirs
# in the two swapped sides above.
rm -f "${RED}/maps/${LONG_CHR}.map" "${RED}/maps/${LONG_CHR}.rfmix.map"
check_eq "setup.sh refits the model on the rotated labels" "$RED_RC" 0
check "and reuses the panel's genotypes (no chromosome read again)" lacks 'reading the panel' "$(cat "${CASE_TMP}/red_setup.log")"
mkdir -p "${G}/HG002R/vcf"
cp "${G}/HG002T/vcf/HG002T.vcf.gz" "${G}/HG002R/vcf/HG002R.vcf.gz"
cp "${G}/HG002T/vcf/HG002T.vcf.gz.tbi" "${G}/HG002R/vcf/HG002R.vcf.gz.tbi"
ANCESTRY_FINE_DIR="$RED" ANCESTRY_FINE_BOOTSTRAP=0 run_step 26b-ancestry-fine.sh HG002R
check_step_exit "26b-ancestry-fine.sh HG002R (HG002T against the rotated labels)"
check "red: without the map of ${LONG_CHR} the red run paints nothing" has 'No chromosome has 100 genotyped SNVs' "$(cat "$STEP_LOG")"
WEA_R=$(frac HG002R continent WEA)
check "red: with rotated labels West Eurasia is no longer first" \
  test "$(awk -F'\t' '$1 == "continent" {print $2; exit}' "${G}/HG002R/ancestry_fine/HG002R_ancestry_fine.tsv" 2>/dev/null)" != WEA
check "red: nor 70% or more (${WEA_R:-none})" test "$(ge "$WEA_R" 0.7 && echo pass)" != pass
rm -rf "$RED"

# --- the Nextflow module ------------------------------------------------------------------
NF_OUT="${G}/nf-ancestry-fine"
SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,vcf,vcf_index\nHG002T,%s,%s\n' "${G}/HG002T/vcf/HG002T.vcf.gz" "${G}/HG002T/vcf/HG002T.vcf.gz.tbi" > "$SHEET"
T0=$(date +%s)
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "${E2E_WORK}/nf-work-ancestry-fine" \
    --input "$SHEET" \
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools ancestry_fine,html_report \
    --ancestry_fine_panel "$PANEL" \
    --outdir "$NF_OUT" \
    --max_cpus 4 --max_memory 14.GB
) 2>&1 | tee "$STEP_LOG"
check_eq "nextflow run with --tools ancestry_fine,html_report exits 0" "${PIPESTATUS[0]}" 0
echo "Nextflow ancestry_fine on HG002T: $(( $(date +%s) - T0 )) s" >> "$E2E_NOTES"
NF_TSV="${NF_OUT}/HG002T/ancestry_fine/HG002T_ancestry_fine.tsv"
check "the module published the fractions table" test -s "$NF_TSV"
# Keyed by level, code and group: a left_out row (fraction ".") has the same
# code at both levels and differs by its group text.
check "its fractions are the bash step's (within 0.01 each)" awk -F'\t' '
  function same(x, y) { if ((x == ".") != (y == ".")) return 0; if (x == ".") return 1; return x - y <= 0.01 && y - x <= 0.01 }
  $1 == "level" {next}
  {k = $1 "|" $2 "|" $3}
  FNR == NR {a[k] = $4; next}
  {b[k] = $4}
  END {
    for (k in a) if (!(k in b) || !same(a[k], b[k])) {print "differs: " k " " a[k] " " b[k]; n++}
    for (k in b) if (!(k in a)) {print "differs: " k " (only in the bash table)"; n++}
    print "differences: " n + 0; exit (n > 0)
  }' "$NF_TSV" "${G}/HG002T/ancestry_fine/HG002T_ancestry_fine.tsv"
# The painting is deterministic (fixed seeds, the same bytes for any thread
# count), so the module's must be the bash step's, row for row.
check "its painting is the bash step's" \
  cmp "${NF_OUT}/HG002T/ancestry_fine/HG002T_painting.bed" "${G}/HG002T/ancestry_fine/HG002T_painting.bed"
check "the HTML report has the fine-grained ancestry card with its painting" \
  bash -c 'grep -qF "Fine-Grained Ancestry" "$1" && grep -qF "data:image/png;base64," "$1"' _ "${NF_OUT}/HG002T/HG002T_report.html"

# The module on TWO's two-sample gVCF: ANCESTRY_FINE_GENOTYPE stops, the
# pipeline ignores the failure of a report-only step and names it in
# failed_tasks.tsv, and nothing is published for TWO.
NF2="${CASE_TMP}/nf-two"
mkdir -p "$NF2"
printf 'sample,vcf,vcf_index,gvcf,gvcf_index\nTWO,%s,%s,%s,%s\n' "${G}/TWO/vcf/TWO.vcf.gz" "${G}/TWO/vcf/TWO.vcf.gz.tbi" \
  "${G}/TWO/vcf/TWO.g.vcf.gz" "${G}/TWO/vcf/TWO.g.vcf.gz.tbi" > "${NF2}/samplesheet.csv"
(
  cd "$NF2" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "${E2E_WORK}/nf-work-ancestry-fine-two" \
    --input "${NF2}/samplesheet.csv" \
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools ancestry_fine \
    --ancestry_fine_panel "$PANEL" \
    --outdir "${NF2}/out" \
    --max_cpus 4 --max_memory 14.GB
) 2>&1 | tee "$STEP_LOG"
echo "+ exit ${PIPESTATUS[0]}"
check "TWO (Nextflow): failed_tasks.tsv names ANCESTRY_FINE_GENOTYPE" has 'ANCESTRY_FINE_GENOTYPE' "$(cat "${NF2}/failed_tasks.tsv" 2>/dev/null)"
check "TWO (Nextflow): the task says why" \
  has 'holds 2 samples; step 26b reads a one-sample file' "$(cat "${E2E_WORK}"/nf-work-ancestry-fine-two/*/*/.command.err 2>/dev/null)"
check "TWO (Nextflow): no fractions published" test ! -e "${NF2}/out/TWO/ancestry_fine/TWO_ancestry_fine.tsv"

finish
