#!/usr/bin/env bash
# Fixed subset of the differential fuzz, plus binary and odd-name controls.
set -euo pipefail
plant() { printf '%b\n' "${2//\~/}" > "$1"; }
plant probe.md '́rs123~ AG\nHG0~02́ *1/~*4\nHTT́ 17/~19\nREPCŃ 17/~19\npatİent *1/~*4\npatıent *1/~*4\n/mnt/u~seŕ\nx.ts.~net́'
plant $'odd\nname.md' '/home/~synthetic/path'
printf '\377\376' >> probe.md
printf '%s\n' 'clean' > clean.md
printf '\0binary' > binary.md
git init -q
git add -- probe.md clean.md binary.md $'odd\nname.md'
