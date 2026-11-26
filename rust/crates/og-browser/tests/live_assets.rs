//! assets/live, its plugins included, holds no outside URL, inline script,
//! on* attribute or eval. Moved here from tests/test_owngenome_live.py
//! (check 7), with its two planted mistakes, which the scan must catch.

use og_browser::{asset_problems, read_assets, repo_root};

#[test]
fn live_assets_load_nothing_from_outside() {
    let texts = read_assets(&repo_root().join("assets/live")).unwrap();
    for f in ["index.html", "live.js", "plugins/charts.js"] {
        assert!(
            texts.contains_key(f),
            "{f} was not read: {:?}",
            texts.keys()
        );
    }
    assert_eq!(asset_problems(&texts), Vec::<(String, String)>::new());

    let mut planted = texts.clone();
    let page = planted.get_mut("index.html").unwrap();
    assert!(page.contains("</body>"));
    *page = page.replace("</body>", "<script>alert(1)</script></body>");
    assert_eq!(
        asset_problems(&planted),
        [("index.html".to_string(), "<script>".to_string())],
        "control: an inline script in index.html"
    );

    let mut planted = texts.clone();
    planted
        .get_mut("live.js")
        .unwrap()
        .push_str("\nfetch('https://cdn.example.com/x.js');\n");
    assert_eq!(
        asset_problems(&planted),
        [(
            "live.js".to_string(),
            "https://cdn.example.com/x.js".to_string()
        )],
        "control: an outside URL in live.js"
    );
}
