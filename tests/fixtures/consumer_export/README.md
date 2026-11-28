# Site table slice for CI

`sites.tsv.gz` is the site table of step 40 (`og consumer
build-sites`, the format `setup.sh --consumer-sites` writes) cut to the e2e
fixture's regions. The e2e case `consumer-export-1` and the all-tools stub
profile (`conf/test_all.config`) read it.

How it was made:

1. dbSNP b157's GRCh38 VCF (`GCF_000001405.40.gz`, the build
   `CONSUMER_DBSNP_BUILD` names) was read over HTTPS with `tabix` for each
   region of the fixture's `regions.bed` (release `fixture-v6`), keeping the
   records whose rsID is in the 23andMe v5 site list of
   `CONSUMER_SITES_COMMIT`: 7,667 rsIDs.
2. The list's rows of those rsIDs and those dbSNP records went through
   build-sites (then `bin/export_consumer.py`, since ported to `og consumer
   build-sites`): 7,623 sites kept; 38 dropped as not an
   SNV in dbSNP, and 6 pseudoautosomal chrX sites whose chrX record lies
   outside the fixture's chrX slices.
3. Since build-sites takes each strand from UCSC's hg19ToHg38 chain, the
   full table built from the pinned inputs changes one site of the slice:
   rs75030631 (chr5) is dropped, because the chain does not place it at
   its dbSNP position and its GRCh37 base is a dbSNP allele on both
   strands. Its row was removed: 7,622 sites kept, 45 dropped. Every other
   row is the one the full table has.
4. build-sites now writes `## table_format: 2 (...)` as the second line,
   and `og consumer export` refuses a table without it. The same line was
   added to the slice; no row changed.

Check of the GRCh37 side: the 7,623 sites kept in step 2 were looked up in
dbSNP b157's GRCh37 VCF (`GCF_000001405.25.gz`) at their GRCh37 position.
Off chrM, 7,337 of 7,338 sites are there with the same reference base (the other,
on chr5, has no GRCh37 record at that position). The chrM sites carry
dbSNP's rCRS position (see docs/40-export-consumer.md).

The slice has rsIDs, positions and reference bases only: no genotype of
anyone.

Sources and terms: the 23andMe v5 site list is from 23andme2vcf
(https://github.com/arrogantrobot/23andme2vcf), Copyright (c) 2015 Rob
Long, MIT licence; its notice is `assets/licences/23andme2vcf-LICENSE.txt`.
dbSNP is NCBI's, with no restriction on use.
