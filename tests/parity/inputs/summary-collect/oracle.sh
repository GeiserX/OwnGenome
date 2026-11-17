#!/bin/sh
# The oracle of og summary collect (tests/parity/oracle.ref): collect() of
# bin/collect_summary.py, given og's options.
#   sh ./oracle.sh PROG --sample S --sample-dir D [--declared-sex X] [--out F] [--assets DIR] [--only GROUP]
# --assets becomes ASSETS_DIR, the global render_report.py --assets sets;
# --only is og's alone (the Python writes every section, and the case
# compares the group's subtrees). The modules collect_summary imports are
# stubbed, except those the group b and c generators stage in pymods/
# and py/. run_delta is stubbed to find no history/ and no recheck/, unless
# py/ brings the real one (the whole-file cases, summary-collect-full.sh).
prog=$1
shift
d=$(mktemp -d)
sed -n 's/^import \([a-z_]*\)  # noqa.*/\1/p' "$prog" | while read -r m; do
  if [ -f "pymods/$m.py" ]; then cp "pymods/$m.py" "$d/"; else : > "$d/$m.py"; fi
done
printf 'def read_recheck(path):\n    return None\n\n\ndef newest_history(sample_dir):\n    return None\n' > "$d/run_delta.py"
# The group c cases bring the modules the oracle imports, as they were at
# its commit, and the two tables it reads beside itself (py/, written by
# tests/parity/gen/summary-collect-c.sh).
if [ -d py ]; then cp py/* "$d"/; fi
cp "$prog" "$d/collect_summary.py"
PYTHONPATH=$d exec python3 -c '
import sys
import collect_summary as c
a = sys.argv[1:]
for opt in ("--assets", "--only"):
    while opt in a:
        i = a.index(opt)
        if opt == "--assets":
            c.ASSETS_DIR = a[i + 1]
        del a[i:i + 2]
sys.exit(c.main(a))
' "$@"
