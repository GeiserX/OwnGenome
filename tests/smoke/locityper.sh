# shellcheck shell=sh
# shellcheck disable=SC2154  # LOCITYPER_PANEL_URL and _TBI_SHA256 come from versions.env, set by image-smoke.sh
# LOCITYPER_IMAGE row: what setup.sh --locityper-data and step 53 run, on the
# fixture's HG002 reads, checked against HG002's own assembly haplotypes.
#
# The HPRC v2.1 panel (LOCITYPER_PANEL_URL) holds HG002's two assembled
# haplotypes, HG002.1 and HG002.2: they are the truth. The row:
#   1. builds a reference of chr6, chr20 and chr22 from the fixture's, and
#      realigns the fixture's reads over HLA class II, CYP2D6 and chr20:10.0-
#      10.5 Mb to it (Locityper needs the BAM's contigs to be the reference's);
#   2. as setup.sh does: jellyfish k-mer counts, the panel's .tbi checked by
#      LOCITYPER_PANEL_TBI_SHA256, the panel read over https 250 kb either
#      side of each locus, locityper target and augment;
#   3. as step 53 does: preproc on chr20:10.0-10.5 Mb (the fixture has no
#      chr17, where Locityper's own background region is), then genotype;
#   4. aligns each called haplotype to HG002.1 and HG002.2 with locityper
#      align and keeps the better pairing: qv.tsv;
#   5. the red control: every haplotype within QV 20 of either truth
#      haplotype taken out of the panel (HG002's own and the near copies other
#      people carry), the database rebuilt from that panel, genotype again:
#      red_qv.tsv, whose calls cannot reach QV 20;
#   6. leaves gt/ and a copy without CYP2D6's success file (gt_cut/) for the
#      row's checks, which run the shipped og locityper report on both in
#      this image: every locus called from the real res.json.gz, and CYP2D6
#      not called once Locityper's per-locus success file is missing.
# The loci are those of assets/locityper_loci.bed (HLA-DRB1, HLA-DQB1,
# CYP2D6); the fixture has no read over FCGR or KIR.
set -e
T=$(nproc)
cd /out

echo "== 1. reference and reads"
samtools faidx /in/ref.fa chr6 chr20 chr22 > mini.fa && samtools faidx mini.fa
samtools view -u -F 0x900 /in/HG002_slice.bam chr6:31500001-33100000 chr20:10000001-10500000 chr22:42000001-42300000 \
  | samtools collate -O -u - collate | samtools fastq -n -1 r1.fq.gz -2 r2.fq.gz -0 /dev/null -s /dev/null -
minimap2 -t "$T" -a -x sr -R '@RG\tID:HG002\tSM:HG002' mini.fa r1.fq.gz r2.fq.gz 2> minimap2.log \
  | samtools sort -@ 2 -o hg.bam - && samtools index hg.bam
rm -f r1.fq.gz r2.fq.gz

echo "== 2. database"
printf 'chr6\t32484347\t32593538\tHLA-DRB1\nchr6\t32652603\t32669077\tHLA-DQB1\nchr22\t42122284\t42152960\tCYP2D6\n' > loci.bed
jellyfish count --canonical --lower-count 2 --out-counter-len 2 --mer-len 25 --threads "$T" --size 400M -o mini.jf mini.fa
wget -q -O panel.tbi "${LOCITYPER_PANEL_URL}.tbi"
sha256sum panel.tbi | cut -d' ' -f1 > panel_tbi_sha256.txt
[ "$(cat panel_tbi_sha256.txt)" = "$LOCITYPER_PANEL_TBI_SHA256" ] || { echo "panel .tbi sha256 $(cat panel_tbi_sha256.txt), want ${LOCITYPER_PANEL_TBI_SHA256}"; exit 1; }
tabix -h "${LOCITYPER_PANEL_URL}##idx##/out/panel.tbi" chr6:32234347-32919077 chr22:41872284-42402960 | bgzip > panel.vcf.gz
tabix -p vcf panel.vcf.gz
locityper target -d db -r mini.fa -j mini.jf -v panel.vcf.gz -g GRCh38 -L loci.bed -@ "$T" > target.log 2>&1 || { tail -n 30 target.log; exit 1; }
locityper augment -d db -@ "$T" > augment.log 2>&1 || { tail -n 30 augment.log; exit 1; }

echo "== 3. preproc and genotype"
locityper preproc -a hg.bam -r mini.fa -j mini.jf -o bg -b chr20:10000001-10500000 -@ "$T" > preproc.log 2>&1 || { tail -n 30 preproc.log; exit 1; }
locityper genotype -a hg.bam -d db -p bg -o gt -@ "$T" > genotype.log 2>&1 || { tail -n 30 genotype.log; exit 1; }
locityper --version > version.txt 2>&1 || true

# The tags of a PAF line, by name: tag["NM"], tag["qv"] (the value after TYPE:).
PAF_TAGS='function tags(  i, j, t) {
  split("", tag)
  for (i = 13; i <= NF; i++) { j = index($i, ":"); t = substr($i, j + 1); tag[substr($i, 1, j - 1)] = substr(t, index(t, ":") + 1) }
}'
# qv GT_DIR OUT: per locus the call and the QV of each called haplotype to
# the truth haplotype it pairs with (the better of the two pairings); 99 for an
# identical sequence. The sequences come from the full panel (db). Locityper
# writes res.json.gz pretty-printed four spaces a level, so its top-level keys
# are the lines indented by exactly four.
qv() {
  : > "$2"
  for locus in HLA-DRB1 HLA-DQB1 CYP2D6; do
    gzip -dc "$1/loci/${locus}/res.json.gz" | awk '/^    "[a-z_]+": / {
        k = $1; gsub(/[":]/, "", k); v = substr($0, index($0, ": ") + 2); sub(/,$/, "", v); gsub(/"/, "", v); r[k] = v }
      END { print r["genotype"] "\t" r["quality"] "\t" (r["weight_dist"] == "" || r["weight_dist"] == "null" ? "nan" : r["weight_dist"]) }' > "$1.${locus}.res"
    gt=$(cut -f1 "$1.${locus}.res")
    case "$gt" in
      *,*,*) echo "${locus}: no diploid genotype in $1 (${gt})"; return 1 ;;
      ?*,?*) ;;
      *) echo "${locus}: no diploid genotype in $1 (${gt})"; return 1 ;;
    esac
    a=${gt%,*} b=${gt#*,}
    # shellcheck disable=SC2046  # the pairs are haplotype names, with no blank
    locityper align -i "db/loci/${locus}/haplotypes.fa.gz" -o "$1.${locus}.paf" -@ 2 -p \
      $(for x in "$a" "$b"; do for y in HG002.1 HG002.2; do [ "$x" = "$y" ] || echo "$x,$y"; done; done | sort -u) > /dev/null 2>&1
    awk -F'\t' -v locus="$locus" -v a="$a" -v b="$b" -v quality="$(cut -f2 "$1.${locus}.res")" \
        -v wd="$(cut -f3 "$1.${locus}.res")" "$PAF_TAGS"'
      function at(x, y) { return ((x, y) in q) ? q[x, y] : 0 }
      function min(u, w) { return u < w ? u : w }
      BEGIN {
        if (quality !~ /^-?[0-9]/) { print locus ": no quality in res.json.gz" > "/dev/stderr"; bad = 1; exit 1 }
        h[1] = "HG002.1"; h[2] = "HG002.2"; x[1] = a; x[2] = b
        for (i = 1; i <= 2; i++) for (j = 1; j <= 2; j++) if (x[i] == h[j]) q[x[i], h[j]] = 99
      }
      /^#/ { next }
      {
        tags()
        if (tag["NM"] == "0") v = 99
        else if ("qv" in tag) v = min(99, tag["qv"] + 0)
        else { print locus ": no qv tag in " FILENAME > "/dev/stderr"; bad = 1; exit 1 }
        q[$1, $6] = v; q[$6, $1] = v
      }
      END {
        if (bad) exit 1
        p1 = at(a, h[1]); p2 = at(b, h[2])
        if (min(at(a, h[2]), at(b, h[1])) > min(p1, p2)) { p1 = at(a, h[2]); p2 = at(b, h[1]) }
        printf "%s\t%s,%s\t%.1f\t%.1f\t%.1f\t%.1f\t%s\n", locus, a, b, p1, p2, min(p1, p2), quality, wd
      }' "$1.${locus}.paf" >> "$2"
    tail -n 1 "$2"
  done
}
echo "== 4. the calls against HG002's haplotypes (locus, genotype, QV 1, QV 2, lower, quality, weighted distance)"
qv gt qv.tsv

echo "== 5. red control: the panel without HG002 and every haplotype within QV 20 of it"
: > red.bed
for l in HLA-DRB1 HLA-DQB1 CYP2D6; do
  # Every panel haplotype against each of HG002's, all aligned: a minimizer
  # divergence cut-off would skip near copies (one base in 250 changes about
  # 6% of the minimizers, so -D 0.05 kept a QV 24 copy of HG002.2's HLA-DQB1).
  gzip -dc "db/loci/${l}/haplotypes.fa.gz" | awk '/^>/ {n = substr($1, 2); if (n !~ /^HG002\./) print n "\tHG002.1\n" n "\tHG002.2"}' > "pairs.${l}.txt"
  locityper align -i "db/loci/${l}/haplotypes.fa.gz" -o "near.${l}.paf" -P "pairs.${l}.txt" -@ "$T" > "near.${l}.log" 2>&1 \
    || { tail -n 20 "near.${l}.log"; exit 1; }
  # HG002's two, and every haplotype identical to either or within QV 20 of
  # it, are left out of red.<locus>.fa.
  { printf 'HG002.1\nHG002.2\n'
    awk -F'\t' "$PAF_TAGS"'/^#/ { next } { tags(); if (tag["NM"] == "0" || tag["qv"] + 0 >= 20) print $1 "\n" $6 }' "near.${l}.paf"
  } > "near.${l}.txt"
  : > "red.${l}.fa"
  gzip -dc "db/loci/${l}/haplotypes.fa.gz" | awk -v l="$l" -v out="red.${l}.fa" '
    FNR == NR { near[$0] = 1; next }
    /^>/ { name = substr($1, 2); n++; if (!(name in near)) kept++ }
    !(name in near) { print >> out }
    END {
      printf "%s\t%d\t%d\n", l, n, n - kept >> "red_removed.tsv"
      print l, n + 0, "haplotypes,", n - kept, "within QV 20 of HG002 removed"
    }' "near.${l}.txt" -
  awk -v l="$l" -v f="/out/red.${l}.fa" -F'\t' '$4 == l {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" f}' loci.bed >> red.bed
done
grep -h '^>' red.*.fa | grep -c '^>HG002\.' > red_hg002.txt || true
locityper target -d dbred -r mini.fa -j mini.jf -L red.bed -@ "$T" > target_red.log 2>&1 || { tail -n 30 target_red.log; exit 1; }
locityper genotype -a hg.bam -d dbred -p bg -o gtred -@ "$T" > genotype_red.log 2>&1 || { tail -n 30 genotype_red.log; exit 1; }
qv gtred red_qv.tsv

echo "== 6. Locityper's own output for og locityper report (run by the row's checks)"
cp -R gt gt_cut && rm gt_cut/loci/CYP2D6/success
rm -rf collate* hg.bam* mini.jf
