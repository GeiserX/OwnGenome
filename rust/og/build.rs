//! Stamps the binary with the git tree hash of the `rust/` it was built
//! from. `rust/og.sh` computes the hash on the host (the build container
//! has no access to the repository's .git) and passes it as OG_TREE_HASH.
//! A plain `cargo build` without it stamps `unknown`, which no check accepts.

fn main() {
    println!("cargo::rerun-if-env-changed=OG_TREE_HASH");
    let hash = std::env::var("OG_TREE_HASH").unwrap_or_default();
    let ok = (40..=64).contains(&hash.len())
        && hash
            .bytes()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase());
    let hash = if ok { hash } else { "unknown".to_owned() };
    println!("cargo::rustc-env=OG_TREE_HASH={hash}");
    println!(
        "cargo::rustc-env=OG_TARGET={}",
        std::env::var("TARGET").expect("cargo sets TARGET")
    );
}
