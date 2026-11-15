//! Clinical links must be live and metadata must agree, without archive GETs.
use og_fresh::http::{Body, Fetch, GitHub, Headers, Request, Want};
use og_fresh::lookups::{Ctx, DATABASES};
use og_fresh::pins::Pins;
use og_fresh::selftest::FakeRegistry;
use og_fresh::{Lookup, failed};
use std::path::Path;
use std::sync::Mutex;

struct Downloads {
    github: GitHub,
    missing: bool,
    bad_length: bool,
    bad_checksum: bool,
    calls: Mutex<Vec<(String, String)>>,
}
impl Fetch for Downloads {
    fn github_api(&self) -> &GitHub { &self.github }
    fn request(&self, r: &Request<'_>) -> Lookup<(Body, Headers)> {
        self.calls.lock().unwrap().push((r.method.into(), r.url.into()));
        if self.missing { return failed("HTTP 404: release does not exist"); }
        if r.url.ends_with("CHECKSUMS") {
            assert_eq!(r.want, Want::Text);
            let body = if self.bad_checksum { "123 5 wrong.tar.gz\n".into() } else {
                "123 5 homo_sapiens_vep_116_GRCh38.tar.gz\n456 6 homo_sapiens_vep_115_GRCh38.tar.gz\n".into()
            };
            Ok((Body::Text(body), Headers::default()))
        } else {
            assert_eq!(r.method, "HEAD", "large archive body was requested");
            Ok((Body::None, Headers::from_pairs(vec![("Content-Length".into(), if self.bad_length { "1" } else { "7664082467" }.into())])))
        }
    }
}
fn check(missing: bool, bad_length: bool, bad_checksum: bool) {
    let fetch = Downloads { github: GitHub { api: "https://invalid".into(), token: None }, missing, bad_length, bad_checksum, calls: Mutex::new(vec![]) };
    let registry = FakeRegistry { tags: |_| panic!("registry is not needed"), digests: |_| panic!("registry is not needed"), asked: Mutex::new(vec![]) };
    let pins = Pins::parse("VEP_CACHE_RELEASE=\"116\"\nPCGR_VEP_CACHE_RELEASE=\"115\"\nPCGR_DATA_BUNDLE=\"20260620\"\nPCGR_DATA_BYTES=\"7664082467\"\n");
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let common = root.join("scripts/lib/common.sh");
    let c = Ctx { fetch: &fetch, registry: &registry, pins: &pins, versions: Path::new("unused"), setup: Path::new("unused"), common: &common, now: "2026-11-15T00:00:00Z".parse().unwrap() };
    for name in ["VEP cache download", "CPSR cache download", "PCGR bundle download"] {
        let (_, lookup) = DATABASES.iter().find(|(n, _)| *n == name).unwrap_or_else(|| panic!("missing liveness row: {name}"));
        let got = lookup(&c);
        let should_fail = missing || (bad_length && name == "PCGR bundle download") || (bad_checksum && name != "PCGR bundle download");
        assert_eq!(got.is_err(), should_fail, "{name}: {got:?}");
    }
    let calls = fetch.calls.lock().unwrap();
    assert!(calls.iter().any(|(_, url)| url.contains("release-115/")));
}
#[test] fn clinical_links_head_only() { check(false, false, false); }
#[test] fn missing_release_fails() { check(true, false, false); }
#[test] fn changed_bundle_length_fails() { check(false, true, false); }
#[test] fn missing_cache_checksum_fails() { check(false, false, true); }
