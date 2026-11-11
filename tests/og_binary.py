"""og_binary.py: where a test or the demo builder finds OwnGenome's og binary.

OG names it when set (Guard's unit-tests job and the docs builds set it to the
one they build). Otherwise the release build of rust/ is used:
$CARGO_TARGET_DIR/release/og, then rust/target/release/og, then on Linux the
static build scripts/setup.sh --og puts in bin/.og. When none is there the
caller stops with the command that builds it, so a missing og fails the test
instead of skipping it.
"""
import os
import platform
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
    for path in paths:
        if os.access(path, os.X_OK):
            return path
    sys.exit("ERROR: og is not built: run (cd rust && cargo build --release --locked -p og), or set OG")
