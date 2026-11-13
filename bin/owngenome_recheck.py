#!/usr/bin/env python3
"""Deprecated: runs `bin/og recheck` with the same arguments. Kept for crontabs until 1.1.0."""
import os, sys
sys.stderr.write("owngenome_recheck.py is deprecated: call bin/og recheck (docs/recheck.md)\n")
og = os.path.join(os.path.dirname(os.path.abspath(__file__)), "og"); os.execv(og, [og, "recheck", *sys.argv[1:]])
