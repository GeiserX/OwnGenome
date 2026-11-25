#!/usr/bin/env bash
# tests/bench/ancestry-engines.sh — step 26b's engines against RFMix and
# ADMIXTURE, on simulated admixed genomes whose truth is known.
#
# Usage: tests/bench/ancestry-engines.sh <panel_dir> <work_dir>
#   panel_dir  a panel installed by scripts/setup.sh --ancestry-fine-panel
#              (samples.tsv, panel/, maps/, model/); it is only read
#   work_dir   where everything is written; a second run keeps every step
#              an earlier run finished and scores again
#
# What it does, with every number fixed so two runs give the same genomes:
#   1. Holds out 19 fixed painting references whose ancestry is one group
#      (4 per continental group, 3 for AMR; the list and its rule are at
#      HELDOUT below), and builds a bench model without them: the panel
#      VCFs without them, each group's frequencies counted again by og
#      ancestry-fine freqs --exclude at the model's SNVs, and a paint_map
#      without them.
#   2. Simulates 5 genomes per scenario (og ancestry-fine simulate) from the
#      held-out people only:
#        S0  unadmixed, one genome per group (AFR, WEA, CSA, EAS, AMR)
#        S1  AFR 0.8, WEA 0.2, T 6
#        S2  WEA 0.5, AMR 0.4, AFR 0.1, T 12
#        S3  WEA 0.5, CSA 0.5, T 20
#        S4  EAS 0.5, WEA 0.5, T 50
#      T is ancestry events per Morgan, about the generations since the mix.
#   3. Runs each engine on the same files: Beagle phases the genome against
#      the bench panel; og ancestry-fine paint and RFMix paint it; og
#      ancestry-fine project and ADMIXTURE -P give the fractions. The
#      commands are step 26b's (scripts/26b-ancestry-fine.sh), with the
#      bench model and no bootstrap.
#   4. Scores each painting against the truth, over every base of the
#      truth's span (a base no row covers counts as wrong):
#        acc    bp accuracy of the copies of each group (phase-free: both
#               copies right, in either order)
#        r2     dosage r2: the squared correlation of the true and the
#               called copies (0, 1 or 2) of each group, over bases and
#               groups, bp-weighted
#        hap    per-haplotype accuracy after the better of the two phase
#               alignments on each chromosome; acc - hap above 0.02 says
#               phase switches of Beagle cost the painting
#        bkpt   phase-free breakpoints (where the copies change): true ones
#               with no call within 1 cM are missed, called ones with no
#               true one within 1 cM are false
#        segs   stretches per haplotype, as the report draws them
#        T      og's fitted T, beside the simulated one
#        |q-t|  the mean and the largest absolute error of the continental
#               fractions against the truth's bp shares
#      and the wall time and peak memory of each engine.
#
# With the default settings and the benchmark panel, the gate checks og
# against the committed five-genome results. Each scenario allows up to
# 0.01 loss in accuracy and dosage r2, and S0 requires accuracy of at
# least 0.98. The gate exits 1 when any row fails and 0 when all rows pass.
# Other settings report scores for both engines without applying floors.
# BENCH_SHIFT=1 can test the gate on the default settings and panel.
#
# Env:
#   THREADS          threads for each tool (default 4)
#   BENCH_CHROMS     chromosomes to simulate, in chromosome order (default chr9 chr15
#                    chr21: about 85,000 panel SNVs and 3.4 Morgans)
#   BENCH_GENOMES    genomes per scenario (default 5)
#   BENCH_SCENARIOS  scenarios to run (default S0 S1 S2 S3 S4)
#   BENCH_SHIFT=1    scores og's painting with each group label moved one
#                    place (AFR read as the next group, and so on): the red
#                    control, which must fail the gate
#   OG               the og Linux binary (default bin/.og/og-linux-<arch>)
# Every tool runs in its image from versions.env, with --cpus THREADS and
# the lowest CPU weight; nothing reads a personal genome. The panel is
# gnomAD's HGDP + 1000 Genomes (CC0), so the genomes and the numbers can be
# published.
set -euo pipefail

PANEL_IN=${1:?Usage: $0 <panel_dir> <work_dir>}
WORK_IN=${2:?Usage: $0 <panel_dir> <work_dir>}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
DEFAULT_THREADS=4
DEFAULT_CHROMS="chr9 chr15 chr21"
DEFAULT_GENOMES=5
DEFAULT_SCENARIOS="S0 S1 S2 S3 S4"
THREADS=${THREADS:-$DEFAULT_THREADS}
CHROMS=${BENCH_CHROMS:-$DEFAULT_CHROMS}
GENOMES=${BENCH_GENOMES:-$DEFAULT_GENOMES}
SCENARIOS=${BENCH_SCENARIOS:-$DEFAULT_SCENARIOS}
SHIFT=${BENCH_SHIFT:-0}
OG=${OG:-${REPO}/bin/.og/og-linux-$(uname -m)}
GROUPS_ALL="AFR WEA CSA EAS AMR"
# HELDOUT: the people the genomes are built from, sample and group. S0
# calls each of its genomes 100% one group, so each person must be that:
# a painting reference whose own-group fraction from og ancestry-fine
# project, against the full panel's frequencies, is at least 0.99, taken
# one population after another. AMR has 3, from Karitiana and Surui only:
# the panel's AMR frequencies are dominated by Maya, so Maya, Pima and
# Colombian people read as about 0.99 AMR against them, while against
# Karitiana and Surui as the AMR group they carry 26 to 56% of other
# groups. A list chosen once keeps the genomes the same however the
# panel's painting map is capped.
HELDOUT="HGDP01405 AFR
HGDP00993 AFR
HGDP00453 AFR
HG02922 AFR
HGDP00995 AMR
HGDP00838 AMR
HGDP01009 AMR
NA20847 CSA
HG03717 CSA
HG03643 CSA
NA20849 CSA
HG00759 EAS
NA18525 EAS
HG00403 EAS
HGDP01309 EAS
HGDP01357 WEA
HGDP01149 WEA
NA06984 WEA
HGDP00511 WEA"
# hold_count <group>: how many of HELDOUT are of that group.
hold_count() { case $1 in AMR) echo 3 ;; *) echo 4 ;; esac; }

# shellcheck source=../../versions.env
. "${REPO}/versions.env"
for v in BCFTOOLS_IMAGE BEAGLE_IMAGE RFMIX_IMAGE ADMIXTURE_IMAGE PLINK2_IMAGE OG_IMAGE; do
  [ -n "${!v:-}" ] || { echo "ERROR: ${v} is not set in versions.env" >&2; exit 2; }
done
case "$THREADS" in ''|*[!0-9]*|0) echo "ERROR: THREADS must be a whole number above 0, got '${THREADS}'" >&2; exit 2 ;; esac
case "$GENOMES" in ''|*[!0-9]*|0) echo "ERROR: BENCH_GENOMES must be a whole number above 0, got '${GENOMES}'" >&2; exit 2 ;; esac
PANEL=$(cd "$PANEL_IN" && pwd)
for f in samples.tsv model/paint_map.tsv model/ref.bim; do
  [ -s "${PANEL}/${f}" ] || { echo "ERROR: ${PANEL}/${f} is missing: install the panel with scripts/setup.sh --ancestry-fine-panel" >&2; exit 2; }
done
for c in $CHROMS; do
  for f in "panel/${c}.vcf.gz" "maps/${c}.map" "maps/${c}.rfmix.map"; do
    [ -s "${PANEL}/${f}" ] || { echo "ERROR: ${PANEL}/${f} is missing" >&2; exit 2; }
  done
done
[ -x "$OG" ] || { echo "ERROR: no og binary at ${OG}: build it with scripts/setup.sh --og, or set OG" >&2; exit 2; }
mkdir -p "$WORK_IN"
WORK=$(cd "$WORK_IN" && pwd)
mkdir -p "${WORK}/stats" "${WORK}/tmp"

# run <stat> <image> <cmd...>: the command in its image, as root, with the
# panel read-only at /panel, the work folder at /work and og at /og; its
# wall seconds and peak memory (the container's cgroup memory.peak) go to
# stats/<stat>.
run() {
  local stat=$1 image=$2 t0 t1 rc=0
  shift 2
  t0=$(date +%s.%N)
  # shellcheck disable=SC2016  # "$@" and $P belong to the inner sh
  docker run --rm --cpus "$THREADS" --cpu-shares 2 --memory 16g -e TMPDIR=/work/tmp -e "P=/work/stats/${stat}.peak" \
    -v "${PANEL}:/panel:ro" -v "${WORK}:/work" -v "${OG}:/og:ro" -w /work "$image" \
    sh -c '"$@"; rc=$?; cat /sys/fs/cgroup/memory.peak > "$P" 2>/dev/null || echo 0 > "$P"; exit $rc' sh "$@" \
    2> >(grep -v "swap limit capabilities" >&2) || rc=$?
  t1=$(date +%s.%N)
  printf '%s\t%s\n' "$(awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.1f", b - a }')" "$(cat "${WORK}/stats/${stat}.peak")" > "${WORK}/stats/${stat}"
  rm -f "${WORK}/stats/${stat}.peak"
  return "$rc"
}

# done_or <marker> : true when an earlier run finished this step.
finished() { [ -e "${WORK}/$1.done" ]; }
mark() { : > "${WORK}/$1.done"; }

echo "Ancestry engines bench: panel ${PANEL}, work ${WORK}, threads ${THREADS}, chromosomes ${CHROMS}"

# --- 1. Held-out people and the bench model ----------------------------------------
if ! finished heldout; then
  echo "$HELDOUT" | tr ' ' '\t' > "${WORK}/heldout.tsv"
  # Each must be a painting reference of its group in this panel: the bench
  # model leaves them out of the paint_map, and RFMix and og must not have
  # painted with them.
  bad=$(awk -F'\t' 'NR == FNR { if (!/^#/ && NF >= 2) g[$1] = $2; next } g[$1] != $2 { print $1 " (" $2 ")" }' \
    "${PANEL}/model/paint_map.tsv" "${WORK}/heldout.tsv" | tr '\n' ' ')
  [ -z "$bad" ] || { echo "ERROR: not a painting reference of that group in ${PANEL}/model/paint_map.tsv: ${bad}" >&2; exit 2; }
  for g in $GROUPS_ALL; do
    k=$(awk -F'\t' -v g="$g" '$2 == g' "${WORK}/heldout.tsv" | wc -l | tr -d ' ')
    [ "$k" -eq "$(hold_count "$g")" ] || { echo "ERROR: ${g} has ${k} people to hold out, not $(hold_count "$g")" >&2; exit 1; }
  done
  cut -f1 "${WORK}/heldout.tsv" > "${WORK}/heldout.txt"
  mark heldout
fi
echo "[1] Held out $(wc -l < "${WORK}/heldout.txt" | tr -d ' ') people (4 per group, 3 for AMR)"

mkdir -p "${WORK}/panel" "${WORK}/model"
for c in $CHROMS; do
  finished "panel/${c}" && continue
  echo "    ${c}: the panel without them"
  # shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
  run "panel.${c}" "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
    'bcftools view --threads "$2" -S ^/work/heldout.txt -Oz -o "/work/panel/$1.vcf.gz.part" "/panel/panel/$1.vcf.gz"
     mv "/work/panel/$1.vcf.gz.part" "/work/panel/$1.vcf.gz"; bcftools index -f -t "/work/panel/$1.vcf.gz"' _ "$c" "$THREADS"
  mark "panel/${c}"
done
if ! finished model; then
  echo "    the bench model: each group's frequencies at the model's SNVs, without the held-out people"
  # The model's SNVs on the bench chromosomes, and the references left.
  awk -F'\t' -v OFS='\t' -v cs=" ${CHROMS} " 'index(cs, " chr" $1 " ") { print "chr" $1, $4 }' \
    "${PANEL}/model/ref.bim" > "${WORK}/model/positions.tsv"
  awk -F'\t' 'NR == FNR { out[$1] = 1; next } FNR > 1 && $5 == 1 && !($1 in out) { print $1 }' \
    "${WORK}/heldout.txt" "${PANEL}/samples.tsv" > "${WORK}/model/refs.txt"
  # shellcheck disable=SC2016  # $1 and $c belong to the inner bash
  run model.vcf "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
    'for c in $1; do bcftools view -S /work/model/refs.txt -T /work/model/positions.tsv -Ob -o "/work/model/$c.bcf" "/panel/panel/$c.vcf.gz"; done
     for c in $1; do echo "/work/model/$c.bcf"; done > /work/model/bcfs.txt
     bcftools concat -f /work/model/bcfs.txt -Oz -o /work/model/adm.vcf.gz' _ "$CHROMS"
  run model.freqs "$OG_IMAGE" /og ancestry-fine freqs --vcf /work/model/adm.vcf.gz \
    --samples /panel/samples.tsv --exclude /work/heldout.txt --out /work/model
  awk -F'\t' -v OFS='\t' 'NR == FNR { out[$1] = 1; next } !/^#/ && NF >= 2 && !($1 in out) { print $1, $2 }' \
    "${WORK}/heldout.txt" "${PANEL}/model/paint_map.tsv" > "${WORK}/model/paint_map.tsv"
  rm -f "${WORK}"/model/chr*.bcf "${WORK}/model/bcfs.txt" "${WORK}/model/adm.vcf.gz"
  mark model
fi
K=$(awk -F'\t' '$1 == "continent" { print $2 }' "${WORK}/model/levels.tsv")
CODES=$(awk -F'\t' '$1 == "continent" { print $3 }' "${WORK}/model/levels.tsv")
echo "    model: $(wc -l < "${WORK}/model/ref.bim" | tr -d ' ') SNVs, continental groups ${CODES}, $(wc -l < "${WORK}/model/paint_map.tsv" | tr -d ' ') painting references"
for g in $GROUPS_ALL; do
  case ",${CODES}," in *",${g},"*) ;; *) echo "ERROR: the model has no group ${g}" >&2; exit 1 ;; esac
done

# Pin the painting references and model markers used for the floors in both
# the input panel and the prepared model, which may come from a cached run.
GATED=0
if [[ "$CHROMS" == "$DEFAULT_CHROMS" && "$GENOMES" == "$DEFAULT_GENOMES" &&
      "$SCENARIOS" == "$DEFAULT_SCENARIOS" && "$THREADS" == "$DEFAULT_THREADS" &&
      ( "$SHIFT" == 0 || "$SHIFT" == 1 ) &&
      $(cksum < "${PANEL}/model/paint_map.tsv") == "2775196785 11436" &&
      $(cksum < "${PANEL}/model/ref.bim") == "1310136243 3254674" &&
      $(cksum < "${WORK}/model/paint_map.tsv") == "166177856 11188" &&
      $(cksum < "${WORK}/model/ref.bim") == "1786942341 309260" ]]; then
  GATED=1
else
  echo "floors not applicable, ungated run (requires default settings and the benchmark panel and model)"
fi

# --- 2. and 3. Genomes and engines ---------------------------------------------------
# scenario <S> <i>: "MU T" of genome i of scenario S.
scenario() {
  local g
  case $1 in
    S0) g=$(echo "$GROUPS_ALL" | awk -v i="$2" '{ print $(1 + (i - 1) % NF) }'); echo "${g}=1 1" ;;
    S1) echo "AFR=0.8,WEA=0.2 6" ;;
    S2) echo "WEA=0.5,AMR=0.4,AFR=0.1 12" ;;
    S3) echo "WEA=0.5,CSA=0.5 20" ;;
    S4) echo "EAS=0.5,WEA=0.5 50" ;;
    *) echo "ERROR: no scenario ${1}" >&2; return 1 ;;
  esac
}

genome() {
  local s=$1 i=$2 n="${1}_g${2}" mu t seed c d="g/${1}_g${2}"
  read -r mu t <<< "$(scenario "$s" "$i")"
  seed=$(( ${s#S} * 100 + i ))
  mkdir -p "${WORK}/${d}/phase" "${WORK}/${d}/paint" "${WORK}/${d}/rfmix" "${WORK}/${d}/frac" "${WORK}/${d}/adm"
  if ! finished "${d}/sim"; then
    echo "    ${n}: simulating (mu ${mu}, T ${t}, seed ${seed})"
    run "${n}.simulate" "$OG_IMAGE" /og ancestry-fine simulate --panel /panel/panel --maps /panel/maps \
      --people /work/heldout.tsv --mu "$mu" --generations "$t" --chroms "${CHROMS// /,}" --seed "$seed" \
      --sample "$n" --out "/work/${d}/sim.vcf.gz" --truth "/work/${d}/truth.bed" > /dev/null
    run "${n}.index" "$BCFTOOLS_IMAGE" bcftools index -f -t "/work/${d}/sim.vcf.gz"
    mark "${d}/sim"
  fi
  if ! finished "${d}/project"; then
    run "${n}.og-project" "$OG_IMAGE" /og ancestry-fine project --model /work/model --vcf "/work/${d}/sim.vcf.gz" \
      --sample "$n" --bootstrap 0 --threads "$THREADS" --out-dir "/work/${d}/frac" > /dev/null
    mark "${d}/project"
  fi
  if ! finished "${d}/admixture"; then
    # ADMIXTURE -P reads the person's genotypes at the model's SNVs, in the
    # .bim's order, and the frequencies as <prefix>.<K>.P.in.
    awk -F'\t' -v OFS='\t' '{ print $1, $2, $2 }' "${WORK}/model/positions.tsv" > "${WORK}/${d}/adm/snvs.bed1"
    run "${n}.plink2" "$PLINK2_IMAGE" plink2 --vcf "/work/${d}/sim.vcf.gz" --extract bed1 "/work/${d}/adm/snvs.bed1" \
      --set-all-var-ids '@:#' --make-bed --threads "$THREADS" --out "/work/${d}/adm/${n}" > /dev/null
    if ! awk -F'\t' 'NR == FNR { sub(/^chr/, "", $1); want[FNR] = $1 ":" $4 ":" $5 ":" $6; n = FNR; next }
                     { sub(/^chr/, "", $1); if ($1 ":" $4 ":" $5 ":" $6 != want[FNR]) exit 1; m = FNR }
                     END { exit m != n }' "${WORK}/model/ref.bim" "${WORK}/${d}/adm/${n}.bim"; then
      echo "ERROR: ${d}/adm/${n}.bim is not the model's ref.bim (SNVs, order or alleles differ)" >&2
      exit 1
    fi
    cp "${WORK}/model/continent.P" "${WORK}/${d}/adm/${n}.${K}.P.in"
    # shellcheck disable=SC2016  # $1 to $3 belong to the inner sh
    run "${n}.admixture" "$ADMIXTURE_IMAGE" sh -c 'cd "$1" && admixture -s 1 -j"$3" -P "$2.bed" '"$K"' > admixture.log 2>&1' \
      _ "/work/${d}/adm" "$n" "$THREADS"
    mark "${d}/admixture"
  fi
  for c in $CHROMS; do
    finished "${d}/phase/${c}" && continue
    # shellcheck disable=SC2016  # $1 to $4 belong to the inner sh
    run "${n}.beagle.${c}" "$BEAGLE_IMAGE" sh -c 'beagle -Xmx8g gt="$1/sim.vcf.gz" chrom="$2" ref="/work/panel/$2.vcf.gz" \
        map="/panel/maps/$2.map" out="$1/phase/$2" impute=false seed=1 nthreads="$3" > "$1/phase/$2.out" 2>&1 \
        || { tail -n 20 "$1/phase/$2.out"; exit 1; }' _ "/work/${d}" "$c" "$THREADS"
    mark "${d}/phase/${c}"
  done
  if ! finished "${d}/paint"; then
    run "${n}.og-paint" "$OG_IMAGE" /og ancestry-fine paint --phased "/work/${d}/phase" --panel /work/panel \
      --maps /panel/maps --model /work/model --fractions "/work/${d}/frac" --sample "$n" --threads "$THREADS" \
      --out-dir "/work/${d}/paint" > "${WORK}/${d}/paint.log"
    [ ! -s "${WORK}/${d}/paint/failed.txt" ] || { echo "ERROR: og did not paint $(cat "${WORK}/${d}/paint/failed.txt")" >&2; exit 1; }
    mark "${d}/paint"
  fi
  for c in $CHROMS; do
    finished "${d}/rfmix/${c}" && continue
    # shellcheck disable=SC2016  # $1 to $3 belong to the inner sh
    run "${n}.rfmix.${c}" "$RFMIX_IMAGE" sh -c 'rfmix -f "$1/phase/$2.vcf.gz" -r "/work/panel/$2.vcf.gz" \
        -m /work/model/paint_map.tsv -g "/panel/maps/$2.rfmix.map" -o "$1/rfmix/$2" --chromosome="$2" \
        --n-threads="$3" --random-seed=1 > "$1/rfmix/$2.log" 2>&1 || { tail -n 20 "$1/rfmix/$2.log"; exit 1; }' \
      _ "/work/${d}" "$c" "$THREADS"
    mark "${d}/rfmix/${c}"
  done
}

echo "[2] Genomes: ${GENOMES} per scenario of ${SCENARIOS}"
for s in $SCENARIOS; do
  scenario "$s" 1 > /dev/null
  for i in $(seq 1 "$GENOMES"); do genome "$s" "$i"; done
done

# --- 4. Scores ------------------------------------------------------------------------
# score <genome dir> <engine dir> <rotate 0|1>: one line of sums (see the header
# of SCORES below) for one painting against the truth.
score() {
  local d=$1 e=$2 rot=$3
  {
    # Truth rows: a start sets the haplotype's group, the end clears it.
    awk -F'\t' -v OFS='\t' '!/^#/ { print $1, $2, 1, "T", $4, $5; print $1, $3, 0, "TE", $4, "-" }' "${d}/truth.bed"
    # Painting rows: row i holds [spos_i - 1, spos_{i+1} - 1), the last
    # one up to its epos, so the stretches touch as the report draws them.
    for c in $CHROMS; do
      awk -F'\t' -v OFS='\t' -v codes="$CODES" -v rot="$rot" '
        BEGIN { k = split(codes, cs, ","); for (j = 1; j <= k; j++) pos[cs[j]] = j }
        /^#Subpopulation/ { sub(/^#Subpopulation order\/codes: /, ""); n = split($0, f, "\t")
          for (j = 1; j <= n; j++) { split(f[j], kv, "="); g = kv[1]; if (rot) g = cs[pos[g] % k + 1]; code[kv[2]] = g }
          next }
        /^#/ { next }
        { print $1, $2 - 1, 1, "C", code[$7], code[$8]; last = $3; chrom = $1 }
        END { if (chrom != "") print chrom, last, 0, "CE", "-", "-" }' "${e}/${c}.msp.tsv"
    done
  } | LC_ALL=C sort -t$'\t' -k1,1 -k2,2n -k3,3n \
    | awk -F'\t' -v OFS='\t' -v codes="$CODES" -v mapdir="${PANEL}/maps" '
      BEGIN { k = split(codes, cs, ",") }
      function cm(c, p,   lo, hi, mid, n, f, l, xs) {
        if (!(c in nmap)) {
          n = 0; f = mapdir "/" c ".map"
          while ((getline l < f) > 0) { split(l, xs, /[ \t]+/); n++; mp[c, n] = xs[4]; mc[c, n] = xs[3] }
          close(f); nmap[c] = n
        }
        n = nmap[c]
        if (n == 0) return 0
        if (p <= mp[c, 1]) return mc[c, 1]
        if (p >= mp[c, n]) return mc[c, n]
        lo = 1; hi = n
        while (hi - lo > 1) { mid = int((lo + hi) / 2); if (mp[c, mid] <= p) lo = mid; else hi = mid }
        if (mp[c, hi] == mp[c, lo]) return mc[c, lo]
        return mc[c, lo] + (mc[c, hi] - mc[c, lo]) * (p - mp[c, lo]) / (mp[c, hi] - mp[c, lo])
      }
      function key(a, b) { return a < b ? a "/" b : b "/" a }
      function flush() {
        if (chrom == "") return
        hap += direct > swap ? direct : swap
        direct = swap = 0
      }
      function account(len,   tk, ck, g, x, y) {
        if (t0 == "-" && t1 == "-") return
        bp += len
        tk = key(t0, t1); ck = key(c0, c1)
        if (c0 != "-" && tk == ck) pf += len
        direct += len * ((t0 == c0) + (t1 == c1)); swap += len * ((t0 == c1) + (t1 == c0))
        for (g = 1; g <= k; g++) {
          x = (t0 == cs[g]) + (t1 == cs[g]); y = (c0 == cs[g]) + (c1 == cs[g])
          w += len; sx += len * x; sy += len * y; sxx += len * x * x; syy += len * y * y; sxy += len * x * y
        }
        if (ltk != "" && tk != ltk) tb[chrom, ++ntb[chrom]] = cm(chrom, at)
        if (lck != "" && ck != lck && c0 != "-" && lc0 != "-") cb[chrom, ++ncb[chrom]] = cm(chrom, at)
        if (lt0 == "" || t0 != lt0) tseg++
        if (lt1 == "" || t1 != lt1) tseg++
        if (c0 != "-" && (lc0 == "" || c0 != lc0)) cseg++
        if (c1 != "-" && (lc1 == "" || c1 != lc1)) cseg++
        ltk = tk; lck = ck; lt0 = t0; lt1 = t1; lc0 = c0; lc1 = c1
      }
      $1 != chrom { flush(); chrom = $1; chroms[++nc] = chrom; t0 = t1 = c0 = c1 = "-"; ltk = lck = lt0 = lt1 = lc0 = lc1 = ""; p = $2 }
      {
        if ($2 > p) { at = p; account($2 - p) }
        p = $2
        if ($4 == "T") { if ($5 == 0) t0 = $6; else t1 = $6 }
        else if ($4 == "TE") { if ($5 == 0) t0 = "-"; else t1 = "-" }
        else if ($4 == "C") { c0 = $5; c1 = $6 }
        else { c0 = c1 = "-" }
      }
      END {
        flush()
        for (j = 1; j <= nc; j++) {
          c = chroms[j]
          for (a = 1; a <= ntb[c]; a++) { hit = 0; for (b = 1; b <= ncb[c]; b++) if ((tb[c, a] - cb[c, b]) ^ 2 <= 1) { hit = 1; break }; if (!hit) missed++; ntrue++ }
          for (b = 1; b <= ncb[c]; b++) { hit = 0; for (a = 1; a <= ntb[c]; a++) if ((tb[c, a] - cb[c, b]) ^ 2 <= 1) { hit = 1; break }; if (!hit) false_++; ncall++ }
        }
        printf "%.0f\t%.0f\t%.0f\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%d\t%d\t%d\t%d\t%d\n",
          bp, pf, hap, w, sx, sy, sxx, syy, sxy, ntrue, missed + 0, ncall + 0, false_ + 0, tseg, cseg
      }'
}

# fractions <genome dir> <Q file>: the mean and the largest |q - truth| of
# the continental fractions, the truth being each group's bp share.
fractions() {
  awk -F'\t' -v codes="$CODES" -v q="$2" '
    BEGIN { k = split(codes, cs, ","); getline l < q; split(l, v, /[ \t]+/) }
    !/^#/ { bp[$5] += $3 - $2; tot += $3 - $2 }
    END { for (g = 1; g <= k; g++) { e = v[g] - bp[cs[g]] / tot; if (e < 0) e = -e; s += e; if (e > m) m = e }
          printf "%.4f\t%.4f\n", s / k, m }' "$1/truth.bed"
}

stat_of() { cat "${WORK}/stats/$1" 2>/dev/null || printf '0\t0\n'; }

# Sums the stats of the per-chromosome runs: seconds added, peak the largest.
stat_sum() {
  local pre=$1 c
  for c in $CHROMS; do stat_of "${pre}.${c}"; done | awk -F'\t' '{ s += $1; if ($2 > m) m = $2 } END { printf "%.1f\t%d\n", s, m }'
}

SCORES="${WORK}/scores.tsv"
{
  printf 'scenario\tgenome\tengine\tbp\tpf\thap\tw\tsx\tsy\tsxx\tsyy\tsxy\ttrue_bkpt\tmissed\tcalled_bkpt\tfalse\ttruth_segs\tsegs\tT_sim\tT_fit\tq_mae\tq_max\tseconds\tpeak_bytes\n'
  for s in $SCENARIOS; do
    for i in $(seq 1 "$GENOMES"); do
      n="${s}_g${i}"; d="${WORK}/g/${n}"
      read -r _ t <<< "$(scenario "$s" "$i")"
      [ "$s" = S0 ] && t=0
      fit=$(awk -F'\t' '$1 == "T" { print $3 }' "${d}/paint/params.tsv")
      printf '%s\t%s\tog\t%s\t%s\t%s\t%s\t%s\n' "$s" "$i" "$(score "$d" "${d}/paint" "$SHIFT")" "$t" "$fit" \
        "$(fractions "$d" "${d}/frac/${n}_continent.${K}.Q")" "$(stat_of "${n}.og-paint")"
      printf '%s\t%s\trfmix\t%s\t%s\t.\t%s\t%s\n' "$s" "$i" "$(score "$d" "${d}/rfmix" 0)" "$t" \
        "$(fractions "$d" "${d}/adm/${n}.${K}.Q")" "$(stat_sum "${n}.rfmix")"
    done
  done
} > "${SCORES}.part"
mv "${SCORES}.part" "$SCORES"

echo ""
[ "$SHIFT" = 1 ] && echo "RED CONTROL: og's group labels are moved one place before scoring."
echo "Per scenario, ${GENOMES} genomes pooled (og = og ancestry-fine paint and project, rf = RFMix and ADMIXTURE -P):"
awk -F'\t' -v OFS='\t' -v gated="$GATED" '
  BEGIN {
    # Equal-reference run, 200-per-group cap, five genomes per scenario.
    # Values match docs/research/ancestry-engine-benchmark.md.
    acc["S0"] = 0.9981; dosage["S0"] = 0.9988
    acc["S1"] = 0.9614; dosage["S1"] = 0.9704
    acc["S2"] = 0.8764; dosage["S2"] = 0.8829
    acc["S3"] = 0.7227; dosage["S3"] = 0.7469
    acc["S4"] = 0.7817; dosage["S4"] = 0.8049
  }
  NR == 1 { next }
  {
    s = $1; e = $3; if (!(s in seen)) { seen[s] = 1; order[++ns] = s }
    for (j = 4; j <= 18; j++) a[s, e, j] += $j
    tsim[s] = $19; if (e == "og") { tfit[s] += $20; nf[s]++ }
    qm[s, e] += $21; if ($22 > qx[s, e]) qx[s, e] = $22; ng[s, e]++
    sec[s, e] += $23; if ($24 > pk[s, e]) pk[s, e] = $24
  }
  function r2(s, e,   w, mx, my, vx, vy, cv) {
    w = a[s, e, 7]; mx = a[s, e, 8] / w; my = a[s, e, 9] / w
    vx = a[s, e, 10] / w - mx * mx; vy = a[s, e, 11] / w - my * my; cv = a[s, e, 12] / w - mx * my
    return (vx > 0 && vy > 0) ? cv * cv / (vx * vy) : (vx == 0 && vy == 0 ? 1 : 0)
  }
  END {
    printf "%-4s %-6s %-6s %-6s %-6s %-6s %-6s %-11s %-11s %-11s %-6s %-6s %-12s %-13s %-13s %s\n",
      "scen", "acc og", "acc rf", "r2 og", "r2 rf", "hap og", "hap rf", "segs t/og/rf", "miss og/rf/t", "false og/rf",
      "T sim", "T fit", "|q-t| og/rf", "sec og/rf", "MiB og/rf", "gate"
    fails = 0
    for (j = 1; j <= ns; j++) {
      s = order[j]
      ao = a[s, "og", 5] / a[s, "og", 4]; ar = a[s, "rfmix", 5] / a[s, "rfmix", 4]
      ho = a[s, "og", 6] / (2 * a[s, "og", 4]); hr = a[s, "rfmix", 6] / (2 * a[s, "rfmix", 4])
      ro = r2(s, "og"); rr = r2(s, "rfmix")
      gate = "UNGATED"
      if (gated) {
        why = ""
        if (ao < acc[s] - 0.01) why = why " acc"
        if (ro < dosage[s] - 0.01) why = why " r2"
        if (s == "S0" && ao < 0.98) why = why " S0-acc<0.98"
        gate = why == "" ? "PASS" : "FAIL:" why
        if (why != "") fails++
      }
      g = ng[s, "og"]
      printf "%-4s %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-11s %-11s %-11s %-6s %-6.1f %-12s %-13s %-13s %s\n",
        s, ao, ar, ro, rr, ho, hr,
        sprintf("%d/%d/%d", a[s, "og", 17] / g, a[s, "og", 18] / g, a[s, "rfmix", 18] / g),
        sprintf("%d/%d", a[s, "og", 14], a[s, "rfmix", 14]) "/" a[s, "og", 13],
        sprintf("%d/%d", a[s, "og", 16], a[s, "rfmix", 16]),
        tsim[s], tfit[s] / nf[s],
        sprintf("%.3f/%.3f", qm[s, "og"] / g, qm[s, "rfmix"] / g),
        sprintf("%.0f/%.0f", sec[s, "og"] / g, sec[s, "rfmix"] / g),
        sprintf("%.0f/%.0f", pk[s, "og"] / 1048576, pk[s, "rfmix"] / 1048576), gate
    }
    print ""
    print "acc, r2 and hap pool the bases of the genomes; segs are stretches per genome (truth/og/RFMix);"
    print "miss: true breakpoints with no call within 1 cM (og/RFMix/all true); false: called ones with no true one"
    print "within 1 cM (og/RFMix); T fit is og'"'"'s fitted T, the mean of the genomes; |q-t| is the mean"
    print "absolute error of the continental fractions (og project / ADMIXTURE -P); sec and MiB are the painting'"'"'s wall time per genome and peak memory."
    if (gated) print "Gate: og acc and r2 >= committed og scores - 0.01 in every scenario, S0 og acc >= 0.98."
    exit fails > 0
  }' "$SCORES" && rc=0 || rc=$?
echo ""
echo "Scores per genome: ${SCORES}"
[ "$GATED" -eq 1 ] || exit "$rc"
if [ "$rc" -ne 0 ]; then
  echo "GATE: FAIL"
  exit 1
fi
echo "GATE: PASS"
