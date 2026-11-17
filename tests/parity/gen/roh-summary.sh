#!/usr/bin/env bash
# roh-summary.sh: writes the bcftools roh outputs of the roh-summary parity
# cases in the current directory, the same every time: runs on chr1 to
# chr22 inside, across and outside the masked bands, runs on chrX, chrY and
# other contigs, the 5 Mb line from both sides, lengths whose megabases
# round on a half, blank-separated fields and lines that are not runs.
# scripts/ci/oracle-diff.sh runs it.
set -euo pipefail
export LC_ALL=C

{
  printf '# This file was produced by: bcftools roh\n# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality (average fwd-bwd phred score)\n'
  printf 'ST\tS\tchr1\t100\t0\t1.0\n'
  printf 'RG\tS\tchrY\t2781480\t47781479\t45000000\t900\t50\n'
  printf 'RG\tS\tchrX\t1\t9000000\t9000000\t900\t50\n'
  printf 'RG\tS\tchrEBV\t1\t6000000\t6000000\t90\t50\n'
  printf 'RG\tS\tchr1_KI270706v1_random\t1\t6000000\t6000000\t90\t50\n'
  printf 'RG\tS\tchr1\t125200001\t143200000\t18000000\t300\t50\n'
  printf 'RG\tS\tchr9\t42300001\t49300000\t7000000\t200\t50\n'
  printf 'RG\tS\tchr1\t120000001\t145000000\t25000000\t500\t50\n'
  printf 'RG\tS\tchr2\t10000001\t17000000\t7000000\t900\t50\n'
  printf 'RG\tS\tchr13\t1\t25000000\t25000000\t900\t50\n'
  printf 'RG\tS\tchr22\t17000001\t40000000\t23000000\t900\t50\n'
  printf 'RG\tS\tchr3\t1\t5000000\t5000000\t900\t50\n'
  printf 'RG\tS\tchr4\t1\t4999999\t4999999\t900\t50\n'
  printf 'RG\tS\tchr5\t1\t5050000\t5050000\t900\t50\n'
  printf 'RG\tS\tchr6\t1\t5150000\t5150000\t900\t50\n'
  printf 'RG\tS\tchr7\t1\t5250000\t5250000\t900\t50\n'
  printf 'RG  S  chr8  1  6000000  6000000  900  50\n'
  printf '  RG\tS\tchr10\t100\t6000099\t6000000\t9\t9\n'
  printf 'RG\tS\tchr23\t1\t9000000\t9000000\t900\t50\n'
  printf 'RG\tS\tchr21\t1\t9000000\n'
  printf 'RG\tS\n'
  printf '\n'
  printf 'RG\tS\tchr11\t50000001\t60000000\t10000000\t9\t9\n'
  printf 'RG\tS\tchr12\t1\t8000000\t8000000.5\t9\t9\n'
} > roh.txt
printf '# empty\n' > empty.txt
