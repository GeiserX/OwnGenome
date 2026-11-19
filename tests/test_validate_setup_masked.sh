#!/usr/bin/env bash
# validate-setup.sh, check_bam_masked_reference: a BAM aligned to the default
# reference is refused when REF_FASTA is the masked one (setup.sh
# --masked-reference), and a BAM aligned to the masked one is refused when
# REF_FASTA is the default, for both aligners. Masking keeps every @SQ name and
# length, so only the aligner's @PG command line tells them apart:
#   minimap2 (step 02, the Nextflow run) names <reference base>.sr.mmi,
#   bwa-mem2 (step 02a) names the FASTA itself.
# The function is read out of validate-setup.sh and run on hand-made headers,
# without Docker. Each verdict is counted from the [OK]/[WARN]/[FAIL] lines,
# and the two refusals are the red controls of the two passes beside them.
#
# The bases, after them: check_reference_content passes the dictionary of
# NCBI's no-ALT analysis set and of its masked copy under their own names,
# and fails the same dictionary with one base of chr21 changed, each set under
# the other's name, and a .dict of another FASTA than the .fai's. The M5 half
# of check_bam_reference fails a BAM whose @SQ M5s are the masked reference's
# against the default one, and when every M5 matches, the masked check above
# trusts them over the @PG file names.
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0

awk '/^check_bam_masked_reference\(\) \{$/ {on = 1} on {print} on && /^}$/ {exit}' \
  "${REPO}/scripts/validate-setup.sh" > "${WORK}/fn.sh"
if grep -q '_maskedGRC' "${WORK}/fn.sh"; then echo "[PASS] check_bam_masked_reference read from validate-setup.sh"
else echo "[FAIL] check_bam_masked_reference not found in validate-setup.sh"; exit 1; fi

STOCK=/genome/reference/GRCh38_no_alt_analysis_set
MASKED=${STOCK}_maskedGRC
sq=$'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr21\tLN:46709983'
mm2() { printf '%s\n@PG\tID:minimap2\tPN:minimap2\tVN:2.28\tCL:minimap2 -t 16 -a -x sr %s.sr.mmi R1.fq.gz R2.fq.gz\n@PG\tID:samtools\tPN:samtools\tPP:minimap2\tCL:samtools sort -o x.bam' "$sq" "$1"; }
bwa() { printf '%s\n@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:2.2.1\tCL:bwa-mem2 mem -t 16 %s.fasta R1.fq.gz R2.fq.gz' "$sq" "$1"; }
none() { printf '%s\n@PG\tID:DRAGEN\tPN:DRAGEN\tCL:dragen -r /ref' "$sq"; }
nocl() { printf '%s\n@PG\tID:minimap2\tPN:minimap2\tVN:2.28' "$sq"; }
# A read group and FASTQ names that hold _maskedGRC, on the default index.
rgname() { printf '%s\n@PG\tID:minimap2\tPN:minimap2\tCL:minimap2 -a -x sr -R @RG\\tID:s_maskedGRC\\tSM:s_maskedGRC %s.sr.mmi s_maskedGRC_R1.fastq.gz s_maskedGRC_R2.fq.gz' "$sq" "$1"; }

# verdict REF HEADER: OK, WARN, FAIL or none, from the function's output.
verdict() {
  (
    # The function sourced below calls these three.
    # shellcheck disable=SC2329
    pass() { echo "OK $1"; }
    # shellcheck disable=SC2329
    warn() { echo "WARN $1"; }
    # shellcheck disable=SC2329
    fail() { echo "FAIL $1"; }
    # shellcheck disable=SC1091
    . "${WORK}/fn.sh"
    REF_FASTA=$1 check_bam_masked_reference "$2"
  ) | awk '$1 == "OK" || $1 == "WARN" || $1 == "FAIL" {s = s (s ? "," : "") $1} END {print (s ? s : "none")}'
}

expect() {  # expect DESC REF HEADER WANT
  local got
  got=$(verdict "$2" "$3")
  if [ "$got" = "$4" ]; then echo "[PASS] $1: $got"
  else echo "[FAIL] $1: got $got, want $4"; FAILS=$((FAILS + 1)); fi
}

expect "masked REF_FASTA, minimap2 BAM on the masked index" "${MASKED}.fasta" "$(mm2 "$MASKED")" OK
expect "masked REF_FASTA, minimap2 BAM on the default index (refused)" "${MASKED}.fasta" "$(mm2 "$STOCK")" FAIL
expect "masked REF_FASTA, bwa-mem2 BAM on the masked FASTA" "${MASKED}.fasta" "$(bwa "$MASKED")" OK
expect "masked REF_FASTA, bwa-mem2 BAM on the default FASTA (refused)" "${MASKED}.fasta" "$(bwa "$STOCK")" FAIL
expect "default REF_FASTA, minimap2 BAM on the masked index (refused)" "${STOCK}.fasta" "$(mm2 "$MASKED")" FAIL
expect "default REF_FASTA, bwa-mem2 BAM on the masked FASTA (refused)" "${STOCK}.fasta" "$(bwa "$MASKED")" FAIL
expect "default REF_FASTA, minimap2 BAM on the default index: nothing to say" "${STOCK}.fasta" "$(mm2 "$STOCK")" none
expect "masked REF_FASTA, no aligner @PG line: a warning" "${MASKED}.fasta" "$(none)" WARN
expect "default REF_FASTA, no aligner @PG line: nothing to say" "${STOCK}.fasta" "$(none)" none
expect "masked REF_FASTA, minimap2 @PG renamed minimap2.1 by a merge" "${MASKED}.fasta" \
  "$(mm2 "$STOCK" | sed 's/ID:minimap2/ID:minimap2.1/; s/PN:minimap2/PN:mm/')" FAIL
expect "default REF_FASTA, masked minimap2 @PG renamed minimap2.1.1 by two merges (refused)" "${STOCK}.fasta" \
  "$(mm2 "$MASKED" | sed 's/ID:minimap2/ID:minimap2.1.1/; s/PN:minimap2/PN:mm/')" FAIL
expect "masked REF_FASTA, an aligner @PG line without CL: a warning" "${MASKED}.fasta" "$(nocl)" WARN
expect "default REF_FASTA, an aligner @PG line without CL: nothing to say" "${STOCK}.fasta" "$(nocl)" none
expect "default REF_FASTA, _maskedGRC only in the read group and FASTQ names: nothing to say" "${STOCK}.fasta" \
  "$(rgname "$STOCK")" none
expect "masked REF_FASTA, _maskedGRC only in the read group and FASTQ names (refused)" "${MASKED}.fasta" \
  "$(rgname "$STOCK")" FAIL
both=$(printf '%s\n%s' "$(mm2 "$STOCK")" "$(bwa "$MASKED" | grep '^@PG')")
expect "masked REF_FASTA, aligner lines naming both references (refused once)" "${MASKED}.fasta" "$both" FAIL
expect "default REF_FASTA, aligner lines naming both references (refused once)" "${STOCK}.fasta" "$both" FAIL

# --- The bases: check_reference_content and the M5 half of check_bam_reference ---
# The dictionary setup.sh makes from NCBI's GCA_000001405.15_GRCh38_no_alt_analysis_set
# (samtools dict 1.24 on the download whose md5 NCBI lists, 2026-11-19): name,
# length and M5 of each of its 195 sequences, in order.
STOCK_DICT=$(cat <<'TSV'
chr1	248956422	6aef897c3d6ff0c78aff06ac189178dd
chr2	242193529	f98db672eb0993dcfdabafe2a882905c
chr3	198295559	76635a41ea913a405ded820447d067b0
chr4	190214555	3210fecf1eb92d5489da4346b3fddc6e
chr5	181538259	a811b3dc9fe66af729dc0dddf7fa4f13
chr6	170805979	5691468a67c7e7a7b5f2a3a683792c29
chr7	159345973	cc044cc2256a1141212660fb07b6171e
chr8	145138636	c67955b5f7815a9a1edfaa15893d3616
chr9	138394717	6c198acf68b5af7b9d676dfdd531b5de
chr10	133797422	c0eeee7acfdaf31b770a509bdaa6e51a
chr11	135086622	1511375dc2dd1b633af8cf439ae90cec
chr12	133275309	96e414eace405d8c27a6d35ba19df56f
chr13	114364328	a5437debe2ef9c9ef8f3ea2874ae1d82
chr14	107043718	e0f0eecc3bcab6178c62b6211565c807
chr15	101991189	f036bd11158407596ca6bf3581454706
chr16	90338345	db2d37c8b7d019caaf2dd64ba3a6f33a
chr17	83257441	f9a0fb01553adb183568e3eb9d8626db
chr18	80373285	11eeaa801f6b0e2e36a1138616b8ee9a
chr19	58617616	85f9f4fc152c58cb7913c06d6b98573a
chr20	64444167	b18e6c531b0bd70e949a7fc20859cb01
chr21	46709983	974dc7aec0b755b19f031418fdedf293
chr22	50818468	ac37ec46683600f808cdd41eac1d55cd
chrX	156040895	2b3a55ff7f58eb308420c8a9b11cac50
chrY	57227415	ce3e31103314a704255f3cd90369ecce
chrM	16569	c68f52674c9fb33aef52dcf399755519
chr1_KI270706v1_random	175055	62def1a794b3e18192863d187af956e6
chr1_KI270707v1_random	32032	78135804eb15220565483b7cdd02f3be
chr1_KI270708v1_random	127682	1e95e047b98ed92148dd84d6c037158c
chr1_KI270709v1_random	66860	4e2db2933ea96aee8dab54af60ecb37d
chr1_KI270710v1_random	40176	9949f776680c6214512ee738ac5da289
chr1_KI270711v1_random	42210	af383f98cf4492c1f1c4e750c26cbb40
chr1_KI270712v1_random	176043	c38a0fecae6a1838a405406f724d6838
chr1_KI270713v1_random	40745	cb78d48cc0adbc58822a1c6fe89e3569
chr1_KI270714v1_random	41717	42f7a452b8b769d051ad738ee9f00631
chr2_KI270715v1_random	161471	b65a8af1d7bbb7f3c77eea85423452bb
chr2_KI270716v1_random	153799	2828e63b8edc5e845bf48e75fbad2926
chr3_GL000221v1_random	155397	3238fb74ea87ae857f9c7508d315babb
chr4_GL000008v2_random	209709	a999388c587908f80406444cebe80ba3
chr5_GL000208v1_random	92689	aa81be49bf3fe63a79bdc6a6f279abf6
chr9_KI270717v1_random	40062	796773a1ee67c988b4de887addbed9e7
chr9_KI270718v1_random	38054	b0c463c8efa8d64442b48e936368dad5
chr9_KI270719v1_random	176845	cd5e932cfc4c74d05bb64e2126873a3a
chr9_KI270720v1_random	39050	8c2683400a4aeeb40abff96652b9b127
chr11_KI270721v1_random	100316	9654b5d3f36845bb9d19a6dbd15d2f22
chr14_GL000009v2_random	201709	862f555045546733591ff7ab15bcecbe
chr14_GL000225v1_random	211173	63945c3e6962f28ffd469719a747e73c
chr14_KI270722v1_random	194050	51f46c9093929e6edc3b4dfd50d803fc
chr14_GL000194v1_random	191469	6ac8f815bf8e845bb3031b73f812c012
chr14_KI270723v1_random	38115	74a4b480675592095fb0c577c515b5df
chr14_KI270724v1_random	39555	c3fcb15dddf45f91ef7d94e2623ce13b
chr14_KI270725v1_random	172810	edc6402e58396b90b8738a5e37bf773d
chr14_KI270726v1_random	43739	fbe54a3197e2b469ccb2f4b161cfbe86
chr15_KI270727v1_random	448248	84fe18a7bf03f3b7fc76cbac8eb583f1
chr16_KI270728v1_random	1872759	369ff74cf36683b3066a2ca929d9c40d
chr17_GL000205v2_random	185591	458e71cd53dd1df4083dc7983a6c82c4
chr17_KI270729v1_random	280839	2756f6ee4f5780acce31e995443508b6
chr17_KI270730v1_random	112551	48f98ede8e28a06d241ab2e946c15e07
chr22_KI270731v1_random	150754	8176d9a20401e8d9f01b7ca8b51d9c08
chr22_KI270732v1_random	41543	d837bab5e416450df6e1038ae6cd0817
chr22_KI270733v1_random	179772	f1fa05d48bb0c1f87237a28b66f0be0b
chr22_KI270734v1_random	165050	1d17410ae2569c758e6dd51616412d32
chr22_KI270735v1_random	42811	eb6b07b73dd9a47252098ed3d9fb78b8
chr22_KI270736v1_random	181920	2ff189f33cfa52f321accddf648c5616
chr22_KI270737v1_random	103838	2ea8bc113a8193d1d700b584b2c5f42a
chr22_KI270738v1_random	99375	854ec525c7b6a79e7268f515b6a9877c
chr22_KI270739v1_random	73985	760fbd73515fedcc9f37737c4a722d6a
chrY_KI270740v1_random	37240	69e42252aead509bf56f1ea6fda91405
chrUn_KI270302v1	2274	ee6dff38036f7d03478c70717643196e
chrUn_KI270304v1	2165	9423c1b46a48aa6331a77ab5c702ac9d
chrUn_KI270303v1	1942	2cb746c78e0faa11e628603a4bc9bd58
chrUn_KI270305v1	1472	f9acea3395b6992cf3418b6689b98cf9
chrUn_KI270322v1	21476	7d459255d1c54369e3b64e719061a5a5
chrUn_KI270320v1	4416	d898b9c5a0118e76481bf5695272959e
chrUn_KI270310v1	1201	af6cb123af7007793bac06485a2a20e9
chrUn_KI270316v1	1444	6adde7a9fe7bd6918f12d0f0924aa8ba
chrUn_KI270315v1	2276	ecc43e822dc011fae1fcfd9981c46e9c
chrUn_KI270312v1	998	26499f2fe4c65621fd8f4ecafbad31d7
chrUn_KI270311v1	12399	59594f9012d8ce21ed5d1119c051a2ba
chrUn_KI270317v1	37690	cd4b1fda800f6ec9ea8001994dbf6499
chrUn_KI270412v1	1179	7bb9612f733fb7f098be79499d46350c
chrUn_KI270411v1	2646	fc240322d91d43c04f349cc58fda3eca
chrUn_KI270414v1	2489	753e02ef3b1c590e0e3376ad2ebb5836
chrUn_KI270419v1	1029	58455e7a788f0dc82034d1fb109f6f5c
chrUn_KI270418v1	2145	1537ec12b9c58d137a2d4cb9db896bbc
chrUn_KI270420v1	2321	bac954a897539c91982a7e3985a49910
chrUn_KI270424v1	2140	747c8f94f34d5de4ad289bc604708210
chrUn_KI270417v1	2043	cd26758fda713f9c96e51d541f50c2d0
chrUn_KI270422v1	1445	3fce80eb4c0554376b591699031feb56
chrUn_KI270423v1	981	bdf5a85c001731dccfb150db2bfe58ac
chrUn_KI270425v1	1884	665a46879bbb48294b0cfa87b61e71f6
chrUn_KI270429v1	1361	ee8962dbef9396884f649e566b78bf06
chrUn_KI270442v1	392061	796289c4cda40e358991f9e672490015
chrUn_KI270466v1	1233	530b7033716a5d72dd544213c513fd12
chrUn_KI270465v1	1774	bb1b2323414425c46531b3c3d22ae00d
chrUn_KI270467v1	3920	db34e0dc109a4afd499b5ec6aaae9754
chrUn_KI270435v1	92983	1655c75415b9c29e143a815f44286d05
chrUn_KI270438v1	112505	e765271939b854dd6826aa764e930c87
chrUn_KI270468v1	4055	0a603090f06108ed7aff75df0767b822
chrUn_KI270510v1	2415	cd7348b3b5d9d0dfef6aed2af75ce920
chrUn_KI270509v1	2318	1cdeb8c823d839e1d1735b5bc9a14856
chrUn_KI270518v1	2186	3fd898b62ca859f50fb8b83e7706352b
chrUn_KI270508v1	1951	7d42a358d472b9cbdfdf30c8742473d0
chrUn_KI270516v1	1300	1cbaaafbbf016906a5bf5886f5a0ecb7
chrUn_KI270512v1	22689	ba1021c82d1230af856f59079e2f71b4
chrUn_KI270519v1	138126	8d754e9c9afd904fba0a2cd577fcc9a1
chrUn_KI270522v1	5674	070b4678e22501029c2e3297115216bc
chrUn_KI270511v1	8127	907ca34a4a2a6673632ebaf513a4c1a4
chrUn_KI270515v1	6361	dd7527ee8e0bdb0a43ca0b2a5456c8c3
chrUn_KI270507v1	5353	311894d0a815eb07c5cac49da851cb4a
chrUn_KI270517v1	3253	913440c38d95c618617ca69bb9296170
chrUn_KI270529v1	1899	4caf890f2586daab8e4b2e2db904f05f
chrUn_KI270528v1	2983	d75c9235f0b8c449fc4352997c56b086
chrUn_KI270530v1	2168	04549369e1197c626669a10164613635
chrUn_KI270539v1	993	19e3a982e67eafef39c5a3e4163f1e17
chrUn_KI270538v1	91309	d60b72221cc7af871f2c757577e4c92a
chrUn_KI270544v1	1202	e62a14b14467cdf48b5a236c66918f0f
chrUn_KI270548v1	1599	866b0db8e9cf66208c2c064bd09ce0a2
chrUn_KI270583v1	1400	b127e2e6dbe358ff192b271b8c6ee690
chrUn_KI270587v1	2969	36be47659719f47b95caa51744aa8f70
chrUn_KI270580v1	1553	1df30dae0f605811d927dcea58e729fc
chrUn_KI270581v1	7046	9f26945f9f9b3f865c9ebe953cbbc1a9
chrUn_KI270579v1	31033	fe62fb1964002717cc1b034630e89b1f
chrUn_KI270589v1	44474	211c215414693fe0a2399cf82e707e03
chrUn_KI270590v1	4685	e8a57f147561b361091791b9010cd28b
chrUn_KI270584v1	4513	d93636c9d54abd013cfc0d4c01334032
chrUn_KI270582v1	6504	6fd9804a7478d2e28160fe9f017689cb
chrUn_KI270588v1	6158	37ffa850e69b342a8f8979bd3ffc77d4
chrUn_KI270593v1	3041	f4a5bfa203e9e81acb640b18fb11e78e
chrUn_KI270591v1	5796	d6af509d69835c9ac25a30086e5a4051
chrUn_KI270330v1	1652	c2c590706a339007b00c59e0b8937e78
chrUn_KI270329v1	1040	f023f927ae84c5cc48dc4dce11ba90f2
chrUn_KI270334v1	1368	53afe12d1371f250a3d1de655345d374
chrUn_KI270333v1	2699	57baf650c47bba9b3a8b7c6d0fb55ad6
chrUn_KI270335v1	1048	eb27188639503b524d2659a23b8262ea
chrUn_KI270338v1	1428	301ef75a6b2996d745eb3464bd352b57
chrUn_KI270340v1	1428	56b462bac20d385cdfcde0155fe4c3a1
chrUn_KI270336v1	1026	69ad2d85d870c8b0269434581e86e30e
chrUn_KI270337v1	1121	16fc8d71a2662a6cfec7bdeec3d810c6
chrUn_KI270363v1	1803	6edd17a912f391022edbc192d49f2489
chrUn_KI270364v1	2855	6ff66a8e589ca27d93b5bac0e5b13a87
chrUn_KI270362v1	3530	bc82401ffd9a5ae711fa0ea34da8d2f0
chrUn_KI270366v1	8320	44a0b65b7ba6bcff37eca202e7d966ea
chrUn_KI270378v1	1048	fc13bda7dbd914c92fb7e49489d1350f
chrUn_KI270379v1	1045	3218bef25946cd95de585dfc7750f63b
chrUn_KI270389v1	1298	2c9b08c57c27e714d4d5259fd91b6983
chrUn_KI270390v1	2387	7a64d89ea14990c16d20f4d6e7283e10
chrUn_KI270387v1	1537	22a12462264340c25e912b8485cdfa91
chrUn_KI270395v1	1143	7c03ca4756c1620f318fb189214388d8
chrUn_KI270396v1	1880	9069bed3c2efe7cc87227d619ad5816f
chrUn_KI270388v1	1216	76f9f3315fa4b831e93c36cd88196480
chrUn_KI270394v1	970	d5171e863a3d8f832f0559235987b1e5
chrUn_KI270386v1	1788	b9b1baaa7abf206f6b70cf31654172db
chrUn_KI270391v1	1484	1fa5cf03b3eac0f1b4a64948fd09de53
chrUn_KI270383v1	1750	694d75683e4a9554bcc1291edbcaee43
chrUn_KI270393v1	1308	3724e1d70677d6b5c4bcf17fd40da111
chrUn_KI270384v1	1658	b06e44ea15d0a57618d6ca7d2e6ac5d2
chrUn_KI270392v1	971	59b3ca8de65fb171683f8a06d3b4bf0d
chrUn_KI270381v1	1930	2a9297cfd3b3807195ab9ad07e775d99
chrUn_KI270385v1	990	112a8b1df94ef0498a0bfe2d2ea5cc23
chrUn_KI270382v1	4215	e7085cdcee6ad62f359744e13d3209fc
chrUn_KI270376v1	1136	59e8fc80b78d62325082334b43dffdba
chrUn_KI270374v1	2656	dbc92c9a92e558946e58b4909ec95dd5
chrUn_KI270372v1	1650	53a9d5e8fd28bce5da5efcfd9114dbf2
chrUn_KI270373v1	1451	b174fe53be245a840cd6324e39b88ced
chrUn_KI270375v1	2378	d678250c97e9b94aa390fa46e70a6d83
chrUn_KI270371v1	2805	a0af3d778dfeb7963e8e6d84c0c54fba
chrUn_KI270448v1	7992	0f40827c265cb813b6e723da6c9b926b
chrUn_KI270521v1	7642	af5bef7cefec7bd7efa729ac6c5be088
chrUn_GL000195v1	182896	5d9ec007868d517e73543b005ba48535
chrUn_GL000219v1	179198	f977edd13bac459cb2ed4a5457dba1b3
chrUn_GL000220v1	161802	fc35de963c57bf7648429e6454f1c9db
chrUn_GL000224v1	179693	d5b2fc04f6b41b212a4198a07f450e20
chrUn_KI270741v1	157432	86eaea8a15a3950e37442eaaa5c9dc92
chrUn_GL000226v1	15008	1c1b2cd1fccbc0a99b6a447fa24d1504
chrUn_GL000213v1	164239	9d424fdcc98866650b58f004080a992a
chrUn_KI270743v1	210658	3b62d9d3100f530d509e4efebd98502c
chrUn_KI270744v1	168472	e90aee46b947ff8c32291a6843fde3f9
chrUn_KI270745v1	41891	1386fe3de6f82956f2124e19353ff9c1
chrUn_KI270746v1	66486	c470486a0a858e14aa21d7866f83cc17
chrUn_KI270747v1	198735	62375d812ece679c9fd2f3d08d4e22a4
chrUn_KI270748v1	93321	4f6c6ab005c852a4352aa33e7cc88ded
chrUn_KI270749v1	158759	c899a7b4e911d371283f3f4058ca08b7
chrUn_KI270750v1	148850	c022ba92f244b7dc54ea90c4eef4d554
chrUn_KI270751v1	150742	1b758bbdee0e9ca882058d916cba9d29
chrUn_KI270752v1	27745	e0880631848337bd58559d9b1519da63
chrUn_KI270753v1	62944	25075fb2a1ecada67c0eb2f1fe0c7ec9
chrUn_KI270754v1	40191	fe9e16233cecbc244f06f3acff3d03b8
chrUn_KI270755v1	36723	4a7da6a658955bd787af8add3ccb5751
chrUn_KI270756v1	79590	2996b120a5a5e15dab6555f0bf92e374
chrUn_KI270757v1	71251	174c73b60b41d8a1ef0fbaa4b3bdf0d3
chrUn_GL000214v1	137718	46c2032c37f2ed899eb41c0473319a69
chrUn_KI270742v1	186739	2f31c013a4a8301deb8ab7ed1ca1cd99
chrUn_GL000216v2	176608	725009a7e3f5b78752b68afa922c090c
chrUn_GL000218v1	161147	1d708b54644c26c7e01c2dad5426d38c
chrEBV	171823	6743bd63b3ff2b5b8985d8933c53290a
TSV
)
# The seven sequences whose M5 the masked copy changes (setup.sh
# --masked-reference: the FASTA with sha256 GRC_MASKED_REF_SHA256, through
# samtools dict 1.24).
MASKED_M5=$'chr9\t1b79085d423b806957b7564497cac5e4\nchr13\t787e7eb2d9187bbc20334062332569d4\nchr16\t9adbaf8ef0094c71470e87eb18e9b5d4\nchr21\t1e771e25108d2a54f1b53661a8d5b6fd\nchrX\t3b4fac479100f43b53358cd77d50085e\nchr11_KI270721v1_random\t894a13477247451089ea90be0365a2f4\nchr22_KI270734v1_random\t7d20b2eb84327e448478329f714f5546'
# chr21 with one base changed (C to T at 30,000,001), through samtools dict 1.24.
CHR21_ONE_BASE=70b7c166708394d0e2601600d67b6b4a

{
  for f in check_reference_content check_bam_reference check_bam_masked_reference; do
    awk -v f="$f" '$0 == f "() {" {on = 1} on {print} on && /^}$/ {exit}' "${REPO}/scripts/validate-setup.sh"
  done
  awk '/^_digest\(\) \{$/ {on = 1} on {print} on && /^}$/ {exit}' "${REPO}/scripts/lib/common.sh"
} > "${WORK}/bases.sh"
for f in check_reference_content check_bam_reference _digest; do
  if grep -q "^${f}() {" "${WORK}/bases.sh"; then echo "[PASS] ${f} read from the scripts"
  else echo "[FAIL] ${f} not found in validate-setup.sh or lib/common.sh"; FAILS=$((FAILS + 1)); fi
done

# with_m5 TABLE CHANGES: TABLE with the M5 of each "name<TAB>M5" line of CHANGES.
with_m5() {
  awk -F'\t' -v OFS='\t' 'NR == FNR { m[$1] = $2; next } ($1 in m) { $3 = m[$1] } { print }' <(printf '%s\n' "$2") <(printf '%s\n' "$1")
}
# dict TABLE: a sequence dictionary as samtools dict writes it.
dict() {
  printf '@HD\tVN:1.0\tSO:unsorted\n'
  awk -F'\t' '{ printf "@SQ\tSN:%s\tLN:%s%s\tUR:file:///genome/ref.fasta\n", $1, $2, ($3 == "" ? "" : "\tM5:" $3) }' <<< "$1"
}

# content NAME TABLE [FAI_TABLE]: verdict of check_reference_content on a
# reference called NAME.fasta whose .dict is TABLE and .fai FAI_TABLE (TABLE
# by default). The FASTA is empty: the check reads the .dict, not the bases.
content() {
  local d="${WORK}/content" ref
  rm -rf "$d"; mkdir -p "$d"
  ref="${d}/$1.fasta"
  : > "$ref"
  awk -F'\t' '{ printf "%s\t%s\t0\t60\t61\n", $1, $2 }' <<< "${3:-$2}" > "${ref}.fai"
  dict "$2" > "${d}/$1.dict"
  touch -r "$ref" "${d}/$1.dict"
  [ -z "${OLD_DICT:-}" ] || touch -t 202001010000 "${d}/$1.dict"
  (
    # shellcheck disable=SC2329
    pass() { echo "OK $1"; }
    # shellcheck disable=SC2329
    warn() { echo "WARN $1"; }
    # shellcheck disable=SC2329
    fail() { echo "FAIL $1"; }
    # shellcheck disable=SC2329
    info() { echo "INFO $1"; }
    # shellcheck disable=SC1091
    . "${WORK}/bases.sh"
    GENOME_DIR=$d REF_FASTA=$ref REF_DICT="${d}/$1.dict" check_reference_content
  ) > "${WORK}/content.out"
  awk '$1 ~ /^(OK|WARN|FAIL|INFO)$/ {s = s (s ? "," : "") $1} END {print (s ? s : "none")}' "${WORK}/content.out"
}

expect_content() {  # expect_content DESC WANT NAME TABLE [FAI_TABLE] [PATTERN]
  local got
  got=$(content "$3" "$4" "${5:-}")
  if [ "$got" = "$2" ] && { [ -z "${6:-}" ] || grep -Eq "$6" "${WORK}/content.out"; }; then echo "[PASS] $1: $got"
  else echo "[FAIL] $1: got $got, want $2${6:+ and a line matching /$6/}"; sed 's/^/       /' "${WORK}/content.out"; FAILS=$((FAILS + 1)); fi
}

DEF=GRCh38_no_alt_analysis_set
MSK=${DEF}_maskedGRC
MASKED_DICT=$(with_m5 "$STOCK_DICT" "$MASKED_M5")
ONE_BASE=$(with_m5 "$STOCK_DICT" $'chr21\t'"$CHR21_ONE_BASE")
expect_content "the default reference: its pinned bases" OK "$DEF" "$STOCK_DICT" "" \
  "OK Reference bases: NCBI's GRCh38 no-ALT analysis set \(the MD5 of each of its 195 sequences"
expect_content "the default reference with one base of chr21 changed (refused)" FAIL "$DEF" "$ONE_BASE" "" \
  "FAIL REF_FASTA has the names and lengths of NCBI's GRCh38 no-ALT analysis set but other bases"
expect_content "upper-case M5 in the .dict reads the same" OK "$DEF" "$(awk -F'\t' -v OFS='\t' '{ $3 = toupper($3); print }' <<< "$STOCK_DICT")"
expect_content "the masked reference: its pinned bases" OK "$MSK" "$MASKED_DICT" "" "OK Reference bases: the masked copy"
expect_content "masked bases under the default name (refused)" FAIL "$DEF" "$MASKED_DICT" "" "its name has no _maskedGRC"
expect_content "default bases under the masked name (refused)" FAIL "$MSK" "$STOCK_DICT" "" "holds the bases of the default one"
expect_content "the masked name with one base of chr21 changed (refused)" FAIL "$MSK" "$(with_m5 "$MASKED_DICT" $'chr21\t'"$CHR21_ONE_BASE")" "" \
  "setup.sh --masked-reference"
expect_content "a .dict of another FASTA than the .fai's (refused)" FAIL "$DEF" "$STOCK_DICT" \
  "$(sed '1s/248956422/248956421/' <<< "$STOCK_DICT")" "does not list the names and lengths"
expect_content "a .dict without M5: a warning" WARN "$DEF" "$(cut -f1,2 <<< "$STOCK_DICT")" "" "has no MD5 \(M5\)"
OLD_DICT=1 expect_content "a .dict older than the FASTA: a warning, then the check" WARN,OK "$DEF" "$STOCK_DICT" "" "is older than REF_FASTA"
expect_content "another reference (a three-sequence slice): nothing pinned" INFO "$DEF" "$(head -n 3 <<< "$STOCK_DICT")" "" "no pinned MD5s"

# The M5 tags of the BAM against REF_DICT. Two sequences of the default
# reference and of the masked one; @PG lines as the aligners write them.
PAIR=$(awk -F'\t' '$1 == "chr20" || $1 == "chr21"' <<< "$STOCK_DICT")
MPAIR=$(with_m5 "$PAIR" "$MASKED_M5")
bam_header() {  # bam_header TABLE PG_REF [no-m5]: SO:coordinate, @SQ lines, a minimap2 @PG line
  printf '@HD\tVN:1.6\tSO:coordinate\n'
  awk -F'\t' -v m="${3:-}" '{ printf "@SQ\tSN:%s\tLN:%s%s\n", $1, $2, (m == "no-m5" ? "" : "\tM5:" $3) }' <<< "$1"
  printf '@PG\tID:minimap2\tPN:minimap2\tVN:2.31\tCL:minimap2 -a -x sr %s.sr.mmi R1.fq.gz R2.fq.gz\n' "$2"
}
# bam NAME TABLE HEADER: verdicts of check_bam_reference then
# check_bam_masked_reference, REF_FASTA NAME.fasta with TABLE as .fai and .dict.
bam() {
  local d="${WORK}/bam" ref
  rm -rf "$d"; mkdir -p "$d"
  ref="${d}/$1.fasta"
  awk -F'\t' '{ printf "%s\t%s\t0\t60\t61\n", $1, $2 }' <<< "$2" > "${ref}.fai"
  dict "$2" > "${d}/$1.dict"
  (
    # shellcheck disable=SC2329
    pass() { echo "OK $1"; }
    # shellcheck disable=SC2329
    warn() { echo "WARN $1"; }
    # shellcheck disable=SC2329
    fail() { echo "FAIL $1"; }
    # shellcheck disable=SC1091
    . "${WORK}/bases.sh"
    export REF_FASTA=$ref REF_DICT="${d}/$1.dict"
    check_bam_reference "$3"
    check_bam_masked_reference "$3"
  ) > "${WORK}/bam.out"
  awk '$1 ~ /^(OK|WARN|FAIL)$/ {s = s (s ? "," : "") $1} END {print (s ? s : "none")}' "${WORK}/bam.out"
}
expect_bam() {  # expect_bam DESC WANT NAME TABLE HEADER [PATTERN]
  local got
  got=$(bam "$3" "$4" "$5")
  if [ "$got" = "$2" ] && { [ -z "${6:-}" ] || grep -Eq "$6" "${WORK}/bam.out"; }; then echo "[PASS] $1: $got"
  else echo "[FAIL] $1: got $got, want $2${6:+ and a line matching /$6/}"; sed 's/^/       /' "${WORK}/bam.out"; FAILS=$((FAILS + 1)); fi
}
expect_bam "default REF_FASTA, a BAM whose M5s are the default's" OK,OK,OK "$DEF" "$PAIR" "$(bam_header "$PAIR" "$STOCK")" \
  "OK BAM's @SQ MD5s \(M5\) match the reference dictionary"
expect_bam "default REF_FASTA, a BAM whose M5s are the masked reference's, @PG naming no _maskedGRC (refused)" OK,FAIL,OK \
  "$DEF" "$PAIR" "$(bam_header "$MPAIR" "$STOCK")" \
  "First difference: sequence 2 \(chr21\) has MD5 1e771e25108d2a54f1b53661a8d5b6fd in the BAM and 974dc7aec0b755b19f031418fdedf293 in the reference dictionary"
expect_bam "masked REF_FASTA, a BAM on the masked bases whose @PG names a FASTA without _maskedGRC: the M5s decide" OK,OK,OK \
  "$MSK" "$MPAIR" "$(bam_header "$MPAIR" /data/renamed)"
expect_bam "masked REF_FASTA, a BAM without M5 on the default index: the @PG rule as before (refused)" OK,OK,FAIL \
  "$MSK" "$MPAIR" "$(bam_header "$MPAIR" "$STOCK" no-m5)" "REF_FASTA is the masked reference, but this BAM was aligned to the default one"

# validate-setup.sh runs the check on the sample BAM's header, right after the @SQ check.
# shellcheck disable=SC2016  # the $ is the script's text, not this shell's
if grep -A1 '^        check_bam_reference "\$BAM_HEADER"$' "${REPO}/scripts/validate-setup.sh" \
   | grep -q '^        check_bam_masked_reference "\$BAM_HEADER"$'; then
  echo "[PASS] validate-setup.sh calls check_bam_masked_reference on the BAM header"
else
  echo "[FAIL] validate-setup.sh does not call check_bam_masked_reference after check_bam_reference"; FAILS=$((FAILS + 1))
fi

echo
if [ "$FAILS" -eq 0 ]; then echo "RESULT: ALL PASS"; else echo "RESULT: ${FAILS} FAILED"; fi
[ "$FAILS" -eq 0 ]
