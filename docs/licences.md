# Licences for downloads

OwnGenome's code uses GPL-3.0-or-later. The data and tools below come from their
publishers and carry their own terms. Review the linked terms before downloading,
installing, or using a resource.

## Annotation data

- AnnotSV's annotation bundle includes third-party data with their own terms. See the [AnnotSV manual](https://www.lbgi.fr/AnnotSV/Documentation/README.AnnotSV_latest.pdf) and [Step 5](05-annotsv.md#annotation-data).
- PCGR/CPSR reference data includes third-party databases with their own terms. The [PCGR bundle report](https://sigven.github.io/pcgrdb/grch38/) lists the sources and licences of the newest bundle release. Setup installs release `20260620`, so that list can differ. See [Step 17](17-cpsr.md#2-pcgr-ref-data-bundle).
- Exomiser's hg38 and phenotype bundles include third-party data with their own terms. See the [Exomiser data source notes](https://exomiser.readthedocs.io/en/15.1.0/advanced_analysis.html#pathogenicitysources) and [Step 48](48-exomiser.md#data).
- The REVEL authors make the scores freely available for non-commercial use. The [Zenodo record](https://doi.org/10.5281/zenodo.7072866) lists ODbL 1.0. Contact the maintainers through the [REVEL download page](https://sites.google.com/site/revelgenomics/downloads) for other uses. See [reference setup](00-reference-setup.md#revel-v13-06-gb) for the manual download.

## Coordinate conversion

UCSC allows free non-commercial use of its hg19ToHg38 chain by independent researchers and nonprofit organisations. Other uses require a commercial licence. See the [UCSC chain download terms](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/). The [chip data guide](chip-data-guide.md#prerequisites)
links to the manual download. Setup fetches the default chain for [Step 43](43-archaic.md#data)
and [Step 44](44-ancient.md#the-panel).

## Optional tools

- Illumina's DRAGEN ORA Decompression Software requires acceptance of its end-user licence. See the [Illumina agreement](https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/dragen-decompression/1000000139522_02_Ora_Decompression_SW_EULA.pdf) and [Step 1](01-ora-to-fastq.md#prerequisites).
- The PacBio Software License Agreement permits Paraphase use only with data generated on a PacBio instrument or provided by PacBio. See the [Paraphase licence](https://github.com/PacificBiosciences/paraphase/blob/v4.1.0/LICENSE.md) and [Step 35b](35-paralogs.md#why-only-pacbio-hifi).
