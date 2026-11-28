# Upgrading OwnGenome

This page says how to move an install to a newer release. It has the steps for every upgrade, then one section for each major release, newest first. [CHANGELOG.md](CHANGELOG.md) lists every change, and each release on [GitHub Releases](https://github.com/GeiserX/OwnGenome/releases) carries its entry.

## Every upgrade

1. Note the release you run now, so you can go back to it: `git -C /path/to/OwnGenome describe --tags`.
2. Read the section below for each major release you skip.
3. Check out the new release:

    ```bash
    cd /path/to/OwnGenome
    git fetch --tags origin
    git checkout v1.0.0   # the release you move to
    ```

4. Name your data folder for this shell: `export GENOME_DIR=/path/to/your/data`. `validate-setup.sh` and `run-all.sh` read it from the environment, so it must be exported.
5. Run setup again with the same folder: `./scripts/setup.sh "$GENOME_DIR"`. It keeps what is already installed and current. It installs what the new release pins: the `og` binaries, images, references and databases.
6. Check each sample before you run it: `./scripts/validate-setup.sh <sample>`.
7. Run the sample as before: `./scripts/run-all.sh <sample> <male|female>`. A rerun reuses each finished step.

To go back, check out the release from step 1. A major release can need more than that: its section says what.

## Upgrading to 1.0

OwnGenome 1.0 is the first stable release. It is one Nextflow pipeline from FASTQ, BAM or VCF to the report. OwnGenome's own code is one Rust binary, `og`. The default reference is a new file, so every BAM from 0.x must be aligned again. From 1.0 on, options, config keys and output file names and formats stay as they are until 2.0. A 1.x release can add a column or a file. A rename or a removal waits for 2.0, and CI fails when a stable output breaks. [Output contract](https://geiserx.github.io/OwnGenome/output-contract/) lists what that covers.

These steps move an install of 0.8.2, or of any earlier 0.x release, to 1.0.

### What changed

- **A new name.** The project is now called OwnGenome. The repository, the docs site and the reports carry the new name. Step 3 points your checkout at it.
- **A new reference.** The pipeline aligns to NCBI's GRCh38 no-ALT analysis set, `reference/GRCh38_no_alt_analysis_set.fasta`. 0.x used the Broad `reference/Homo_sapiens_assembly38.fasta`, which has ALT, HLA and decoy contigs. Reads at CYP2D6, the HLA genes and KIR got mapping quality 0 on it, so depth there was wrong. `validate-setup.sh` refuses a BAM aligned to another reference. [Realigning after a reference change](https://geiserx.github.io/OwnGenome/realignment/) says why and what to redo.
- **`run-all.sh` runs Nextflow.** It takes the same arguments as in 0.x, `<sample> <male|female>`. It writes a samplesheet and starts `nextflow run main.nf -resume`, so a full run needs Java 17 or later and Nextflow. The results go to the same `<genome_dir>/<sample>/` folders. Options after the sex go to Nextflow as they are. The bash scheduler is gone. The numbered scripts still run one step each.
- **The `og` binary.** OwnGenome's own tools are one static Rust binary, `bin/og`. Setup installs it. On a release tag it downloads the release's binaries and checks them against `SHA256SUMS`. On `main` it builds them in a container, which takes longer. `./scripts/setup.sh --og` installs only `og`. A step that needs `og` stops and says so when it is missing.
- **No Python of OwnGenome's own.** 0.x ran its own Python inside the Python image for several steps. That code is now part of `og`. The host needs no `python3`.
- **Setup installs more by default.** `./scripts/setup.sh <genome_dir>` now also installs the clinical data. Steps 13, 23, 30 and 31 read the VEP 116 cache. Step 17 reads the VEP 115 cache and the PCGR bundle. 0.x left these to you. `--no-clinical-data` skips them. Each optional step has its own setup flag, and `--everything` installs all of them. [Reference data setup](https://geiserx.github.io/OwnGenome/00-reference-setup/) lists them.
- **The downloads folder.** Setup keeps every raw download in `<genome_dir>/downloads`, or in the folder that `OG_DOWNLOADS` names. A rerun then downloads nothing it already has. The kept files take about 69 GB for a plain setup. `--no-keep-downloads` deletes each one after its install. Setup moves an archive that 0.x left in `vep_cache/` or `pcgr_data/` into the downloads folder, and checks it there.
- **CPSR 2.3.2.** Step 17 runs PCGR/CPSR 2.3.2 from `ghcr.io`, not 2.2.5. It reads the bundle `pcgr_data/20260620/` and the VEP 115 cache. 0.x read the bundle `20250314` and the VEP 113 cache. The new bundle and cache are about 31 GB to download (7 GB and 24 GB). Its variant table is now `cpsr/<sample>.cpsr.grch38.classification.tsv.gz`, not `<sample>.cpsr.grch38.snvs_indels.tiers.tsv`.
- **The Nextflow pipeline starts from reads.** The samplesheet takes `fastq_1`, `fastq_2` and `sex` columns. A 0.x samplesheet (`sample,vcf,vcf_index,bam,bam_index`) still works. An option or a `--tools` name the pipeline does not know stops the run with its name, so a misspelt option is no longer ignored. `--genome` is gone. `cyrius` is no longer in the default `--tools`, because its licence is non-commercial. `run-all.sh` runs it when `setup.sh --cyrius` installed it. `--slivar_bin` is ignored with a warning in 1.0.

### Steps

1. Note the release you run now: `git -C /path/to/OwnGenome describe --tags`. That release is your way back.
2. Install Java 17 or later and the Nextflow release that `NEXTFLOW_VERSION` in `versions.env` pins. [Install the requirements](https://geiserx.github.io/OwnGenome/getting-started/#install-the-requirements) has the commands for each platform, and a check.
3. Point the checkout at the new name, then check out 1.0. GitHub still forwards the old name, but the new one is safer:

    ```bash
    cd /path/to/OwnGenome
    git remote set-url origin https://github.com/GeiserX/OwnGenome.git
    git fetch --tags origin
    git checkout v1.0.0
    ```

4. Name your data folder for this shell: `export GENOME_DIR=/path/to/your/data`. Setup, `validate-setup.sh` and `run-all.sh` all use it.
5. Run setup: `./scripts/setup.sh "$GENOME_DIR"`. It downloads the new reference, installs `og`, pulls the new images and installs the clinical data. It checks the free disk before the first download and stops when the downloads do not fit. Add the flag of each optional step you use, or `--no-clinical-data` if you do not run steps 13, 17, 23, 30 and 31.
6. Align each sample again. [Realigning after a reference change](https://geiserx.github.io/OwnGenome/realignment/) has the steps, one sample at a time, and the check that a BAM is on the new reference. Results made from a 0.x BAM have to be made again.
7. Check each sample: `./scripts/validate-setup.sh <sample>`. Fix each `[FAIL]` line before you run.
8. Run each sample: `./scripts/run-all.sh <sample> <male|female>`. A full run from FASTQ can take more than a day on 8 CPUs. Start it inside `tmux` or `screen`.
9. If you start Nextflow yourself, remove `--genome` from your command, and add `cyrius` to `--tools` if you want it.
10. Once each sample has run on 1.0, delete the files 1.0 no longer reads. They are the old reference, `reference/Homo_sapiens_assembly38.fasta`, with its `.fai`, `.dict`, `.sr.mmi` and any BWA index beside it. They are also `pcgr_data/20250314/` and `vep_cache/homo_sapiens/113_GRCh38/`.

### If you ran `main` before 1.0

Some parts of 1.0 replaced work that only `main` had between 0.8.2 and 1.0. If you ran a checkout of `main` in that time, these steps apply too.

**The database re-check.** `bin/owngenome_recheck.py` is gone. Point each crontab line at `bin/og recheck` with the same arguments. [Database re-check](https://geiserx.github.io/OwnGenome/recheck/) has the lines.

**The live site.** Up to 1.0, `main` served the live site with the Python server `bin/owngenome_live.py`. 1.0 serves it with `og live serve`, which reads the same state folder, so the runs, their messages and the password hash carry over. Switch between runs: a run that is still going keeps its old `live.config`, which has no token.

1. Note the commit you run now. It still has the Python server, and it is your way back: `git -C /path/to/OwnGenome rev-parse HEAD`.
2. Stop the site (`systemctl stop owngenome-live`) and copy what it uses: `cp -a <state dir> <state dir>.before-og` and `cp /etc/systemd/system/owngenome-live.service ~/owngenome-live.service.before-og`. The state folder is `OWNGENOME_LIVE_STATE_DIR`, else `GENOME_DIR/.owngenome-live`. Both servers must use the same one. If yours is another folder, set `Environment=OWNGENOME_LIVE_STATE_DIR=...` in the unit.
3. Check out the release, as in [Every upgrade](#every-upgrade) step 3: `git fetch --tags origin`, then `git checkout v1.0.0`. Then install `og`: `scripts/setup.sh --og`. On the release tag it downloads the release's binaries and checks them.
4. Give the state folder to the account in the unit's `User=`, and let no other account write to it: `chown -R genome <state dir>` and `chmod -R go-w <state dir>`. `og live serve` does not start on a state folder that another account owns or may write to.
5. Remove `--plugin-dir` from `ExecStart` and `plugin_dir=` from the `--config` file. The plugins are part of `og` now, and `og live serve` refuses both. A plugin of your own (`owngenome_live_<name>.py`) does not run under `og`. Port it to `rust/crates/og-live/src/plugins.rs` ([Plugins](https://geiserx.github.io/OwnGenome/live-site/#plugins)), or stay on the Python server until you do. The Run health card and `/metrics` are built in, with the same flags ([Run Alerts](https://geiserx.github.io/OwnGenome/alerts/)).
6. If `ExecStart` runs `bin/owngenome_live.py` itself, point it at `scripts/owngenome-live.sh` with the same flags. Then run `systemctl daemon-reload` and `systemctl restart owngenome-live`.
7. Check the site. `systemctl status owngenome-live` shows `og live serve` running. `<state dir>/ingest-token` exists with mode 600. The site loads in a browser, and the home page lists the earlier runs. Start the next run with `run-all.sh` as before: it sends the token itself.

To go back, stop the site. Check out the commit from step 1 (`git -C /path/to/OwnGenome checkout <commit>`, or `git worktree add` it in another folder). Put back the unit you copied in step 2, run `systemctl daemon-reload` and start the site. The unit's `scripts/owngenome-live.sh` is then that commit's script, which starts the Python server. The Python server reads the same `runs/<run>/events.jsonl` files, including the messages that `og live serve` took. It ignores the files that `og live serve` added (`ingest-token` and the `.lock` files). A 1.0 `run-all.sh` posts to `/ingest/<run>` with the token in a header, and the Python server takes those posts too, so a run that started meanwhile keeps reporting. The copy from step 2 holds the state as it was before the switch. Put it back only if you want that state, because the messages taken since then are only in the current folder.
