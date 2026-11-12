"""og_binary.py: where a test or the demo builder finds OwnGenome's og binary.

OG names it when set (Guard's unit-tests job and the docs builds set it to the
one they build). Otherwise the release build of rust/ is used:
$CARGO_TARGET_DIR/release/og, then rust/target/release/og, then on Linux the
static build scripts/setup.sh --og puts in bin/.og. In a git checkout only a
candidate stamped with this checkout's rust/ tree hash (rust/og.sh tree-hash)
is taken, so a stale build never shadows a fresh one. When none fits, the
caller stops with the command that builds it, so a missing or stale og fails
the test instead of skipping it or running old code.
"""
import os
import platform
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def og_binary():
    """The path of an og binary that runs on this machine."""
    named = os.environ.get("OG")
    if named:
        if not os.access(named, os.X_OK):
            sys.exit(f"ERROR: OG={named} is not an executable og binary")
        return named
    targets = [os.environ.get("CARGO_TARGET_DIR"), os.path.join(REPO, "rust", "target")]
    paths = [os.path.join(t, "release", "og") for t in targets if t]
    if platform.system() == "Linux":
        arch = {"amd64": "x86_64", "arm64": "aarch64"}.get(platform.machine(), platform.machine())
        paths.append(os.path.join(REPO, "bin", ".og", f"og-linux-{arch}"))
    expected = _tree_hash()
    stale = []
    for path in paths:
        if not os.access(path, os.X_OK):
            continue
        if expected is None or _stamp(path) == expected:
            return path
        stale.append(path)
    build = "(cd rust && OG_TREE_HASH=$(./og.sh tree-hash) cargo build --release --locked -p og)"
    if stale:
        sys.exit(f"ERROR: og is stale ({', '.join(stale)} not built from this rust/): run {build}, or set OG")
    sys.exit(f"ERROR: og is not built: run {build}, or set OG")


def _tree_hash():
    """This checkout's rust/ tree hash, or None outside a git checkout."""
    r = subprocess.run(["bash", os.path.join(REPO, "rust", "og.sh"), "tree-hash"],
                       capture_output=True, text=True, check=False)
    return r.stdout.strip() if r.returncode == 0 else None


def _stamp(path):
    """The tree hash stamped into an og binary, or None (a plain cargo build)."""
    with open(path, "rb") as f:
        m = re.search(rb"og-tree-hash:([0-9a-f]{40,64});", f.read())
    return m.group(1).decode() if m else None
