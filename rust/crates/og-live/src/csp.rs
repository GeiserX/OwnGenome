//! The live site's Content-Security-Policy and the headers that go with it
//! on every answer of `og live serve`.
//!
//! [`CSP`] is the policy of the site's own pages: scripts, styles, data and
//! frames from this site only, nothing inline, no outside URL, and no page
//! of another site may frame it. A run's own HTML (the report, CPSR,
//! PharmCAT, MultiQC) carries inline styles and, for third-party pages,
//! scripts and outside links; it is sent under [`RUN_HTML_CSP`], sandboxed:
//! no script runs, nothing is fetched from outside, and it keeps its
//! styles. `allow-same-origin` lets its links reach the run's other files;
//! without `allow-scripts` that grants the page nothing else.

/// The policy of every page, script, style and data answer of the site.
pub const CSP: &str = "default-src 'none'; script-src 'self'; connect-src 'self'; style-src 'self'; \
     img-src 'self' data:; form-action 'self'; base-uri 'none'; frame-src 'self'; frame-ancestors 'none'";

/// The policy of a run's own HTML and SVG files.
pub const RUN_HTML_CSP: &str = "sandbox allow-same-origin allow-popups allow-popups-to-escape-sandbox; \
     default-src 'none'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; form-action 'none'; \
     base-uri 'none'; frame-ancestors 'none'";

/// `X-Frame-Options` of every answer but the genome viewer's.
pub const DENY: &str = "DENY";

/// `X-Frame-Options` of the genome viewer, which the map's reads panel
/// frames.
pub const SAME_ORIGIN: &str = "SAMEORIGIN";

/// The headers every answer carries besides its type and length: the policy
/// `csp`, no type sniffing, no referrer, no cache, and the frame option.
pub fn common(csp: &str, frame_options: &str) -> [(&'static str, String); 5] {
    [
        ("Content-Security-Policy", csp.to_owned()),
        ("X-Content-Type-Options", "nosniff".to_owned()),
        ("Referrer-Policy", "no-referrer".to_owned()),
        ("Cache-Control", "no-store".to_owned()),
        ("X-Frame-Options", frame_options.to_owned()),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_site_policy_allows_nothing_from_outside_and_no_framing() {
        for part in [
            "default-src 'none'",
            "script-src 'self'",
            "frame-ancestors 'none'",
            "base-uri 'none'",
        ] {
            assert!(CSP.contains(part), "{part}");
        }
        assert!(!CSP.contains("unsafe"));
        assert!(RUN_HTML_CSP.starts_with("sandbox ") && !RUN_HTML_CSP.contains("allow-scripts"));
        let h = common(CSP, DENY);
        assert_eq!(h[0], ("Content-Security-Policy", CSP.to_owned()));
        assert_eq!(h[4].1, "DENY");
    }
}
