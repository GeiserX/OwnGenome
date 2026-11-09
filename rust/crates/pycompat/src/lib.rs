//! Python-exact primitives for the Rust ports of OwnGenome's Python tools.
//!
//! Every function here returns the bytes CPython 3.15 returns for the same
//! input. The differential test (`tests/rust/pycompat-diff.sh`) checks that
//! against `python3 -I` in the pinned `PYTHON_IMAGE`.
//!
//! - [`fmt`]: `repr(float)`, the `f`, `e`, `g` and `%` format types, `{:,}`,
//!   zero padding, `ljust`/`rjust`, `round(x)` and `round(x, n)`.
//! - [`num`]: `float(str)`, `int(str)`, `math.fsum` and the builtin `sum`.
//! - [`json`]: `json.dumps` and a `json.loads` with Python's leniency.

pub mod fmt;
pub mod json;
pub mod num;

#[cfg(all(test, feature = "planted-failure"))]
mod planted {
    /// Always fails. CI runs it with `--features planted-failure` and
    /// requires a red result, so a test runner that reports nothing cannot
    /// pass for a green one.
    #[test]
    fn planted_failure() {
        assert_eq!(crate::fmt::repr(0.1), "0.2", "planted failure");
    }
}
