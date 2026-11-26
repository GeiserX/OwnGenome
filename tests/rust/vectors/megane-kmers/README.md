# MEGAnE k-mer set vectors

The oracle of `rust/crates/og-setup/tests/megane_kmers.rs`: small synthetic FASTAs and the `.mk` and `.mi` files that MEGAnE's own `build_kmerset` wrote for them. `og setup megane-kmers` must write the same bytes for each FASTA, on any thread count and with any pass size.

MEGAnE is MIT-licensed (RIKEN). The pinned `MEGANE_IMAGE` in `versions.env` holds MEGAnE at commit `6bfbcff`. The image tag says v1.0.1.beta, but its k-mer builder is the one at that commit. Its output is the same as the tag's, because the only change there is how much memory the builder reserves.

## Files

| File | What it is |
|---|---|
| `rules.fa` | One contig per input rule, each with the copies that make its k-mers occur twice, so the rule shows in the output (list below). The last line has no line end. |
| `random.fa` | Four random chromosomes (lines of 70 and 80 bases) with planted repeats: copies that change one base in 40, reverse-complement copies, CA repeats, poly-A, N runs, a soft-masked stretch and other letters. |
| `twice.fa` | Two contigs, `x` and `y`. |
| `*.fa.fai` | `samtools faidx` of each FASTA. `twice.fa.fai` was edited by hand to list `y` twice, before and after `x`, so every k-mer of `y` counts twice. |
| `*.mk`, `*.mi` | What MEGAnE's `build_kmerset` wrote. |
| `record.sh` | The recording steps. |

The rules `rules.fa` covers:

- Both strands. Each window adds its forward value and the value of its reverse complement (`fwd` and `rev`).
- A palindromic window adds one value. `pal` holds one palindrome twice, so it is kept, and another once, so it is not.
- The coding is A=0 T=1 G=2 C=3, with the first base in the top bits.
- Lowercase reads as uppercase (`lower` and `upper`).
- N and n drop each window that holds them (`nblock` and `nfree`).
- The N-at-base-0 quirk. When the only N in a contig's first 32 bases is base 0, MEGAnE keeps that first window. The N reads as A forward and as T in the reverse complement (`n0`, `a0`, `t0`, `n0lower`). An N at base 0 and at base 5 drops the window (`n0n5`, `a0a5`).
- IUPAC letters and `-` do not drop windows. They read as A forward and as T in the reverse complement (`iupac`, `iupac_a`, `iupac_t`).
- No window spans two contigs (`edge1`, `edge2`, `whole`).
- Contigs of 31, 32 and 33 bases (`s31`, `s32`, `s33`, `shorts`).
- Lines of 7, 50, 60 and 61 bases.
- A k-mer is kept when it occurs 2 or more times.

## Recording

The FASTAs are what the test's generator writes, and the test checks this. To write them again:

```bash
cd rust
OG_WRITE_MEGANE_VECTORS=1 cargo test -p og-setup --test megane_kmers the_committed_fastas
```

Then index them and run MEGAnE's builder. This needs Docker. MEGAnE's image is amd64 only, and emulation is fine for inputs this small:

```bash
tests/rust/vectors/megane-kmers/record.sh
```

`record.sh` runs `samtools faidx` in `SAMTOOLS_IMAGE`, edits `twice.fa.fai`, and runs `build_kmerset -fa <name>.fa -prefix <name> -outdir out` in `MEGANE_IMAGE` for each FASTA. Then it prints the sha256 of every file. The repo's `.gitignore` ignores `*.fa` and `*.fa.fai`, so stage a new recording with `git add -f`, or the FASTAs and their indexes stay out of the commit without a word. The folder must be visible to Docker. On a Docker VM that shares only the home folder, copy `versions.env` and this folder (at the same relative path) under home, and run it there.

The recording of 2026-11-26 gave:

```
3617fb6ec47464ac2f6dbe3173b6f500805bf678392251849ec777bbae5535bb  rules.fa
f254e1a4da6ada93e96d92be391ac03c1db520ab4be12319a6f51e2e50cb0764  random.fa
5bf2755c5e515e89e63130c98e2a3a4d9514651d572a4a2301afae4bc02c002f  twice.fa
d3bccdec049ad603ee12fb3d5cfb3b8d94ba475e5478e7d55067e74c7ac6a325  rules.fa.fai
5c74781fb9510bc1a02a537022d5bf2f78dcc485da0924798adc6536dcda8580  random.fa.fai
f5aead14297ac5be7c72d6b9d502ad59c2d1a312ebc19920edda8214b3e1ad0d  twice.fa.fai
e116e54900e6b32d0c74ede5d3469e0785ecc53909fd5eaa86db8abb693cb4b9  rules.mk    (1392 k-mers)
30df69906d9a76e11ee93bbf2b67c78398e9a15ce67d06eb03826ff1c5c313da  random.mk   (2948 k-mers)
1760c8b80acf1a21f81d13be66ed66582e9e37a346f58a2f76aaae9d27deaa62  twice.mk    (98 k-mers)
a1f34d7e729c9cd1f9fe6494999ee72f206ff6d00d72180b0979784a0e49fdfc  rules.mi
65f68e37473fc9b8ce17d96bc012459aa6b6c04167c1c61367deea5b909c7693  random.mi
4960a9ce196ee53a9ae6f2b038ebf5bb7949312406eabe79835e766b4a3d0d88  twice.mi
```

## Input that og refuses

MEGAnE's builder reads some input in an undefined way, or reads a FASTA that its `.fai` does not describe. og refuses this input, so it has no vector here. The test builds each case itself:

- Windows (CRLF) line ends. MEGAnE reads each `\r` as a base.
- Bytes of 0x80 or more. MEGAnE indexes its table with a negative number.
- A line end inside a contig's bases, or a FASTA shorter than the `.fai` says. Both mean a stale `.fai`.
- A gzip-compressed FASTA.
- A `.fai` line that is not five fields, with whole numbers in columns 2 to 5 and at least 2 bytes per line.
