//! The public demo through Chrome at a Pages-style subpath.

#[path = "live_demo_perf.rs"]
mod perf;

use og_browser::{Browser, Scratch, StaticServer, chrome_or_skip, js_str, repo_root};
use pycompat::json::{Value, loads};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Duration;

const MAIN: &str = "DEMO-001";
const PARTNER: &str = "DEMO-002";
const CHILD: &str = "DEMO-003";
const FAMILY: &str = "DEMO-FAMILY";
const READY: &str = "document.documentElement.dataset.pgpStatus === 'ready'";
const RESULTS: &str = "document.documentElement.dataset.pgpResults";
const MAP: &str = "document.getElementById('genome-map')";
const WAIT: Duration = Duration::from_secs(12);

fn og() -> Option<PathBuf> {
    let path = std::env::var_os("OG").map(PathBuf::from).or_else(|| {
        let test = std::env::current_exe().ok()?;
        Some(test.parent()?.parent()?.join("og"))
    });
    if path.as_ref().is_some_and(|p| p.is_file()) {
        return path;
    }
    assert_ne!(
        std::env::var("CI").as_deref(),
        Ok("true"),
        "build og or set OG in CI"
    );
    eprintln!("[SKIP] build og or set OG to check the public demo");
    None
}

fn build(og: &Path, site: &Path, work: &Path) {
    let built = Command::new(og)
        .args(["live", "demo", "--records", "2000", "--out"])
        .arg(site)
        .arg("--work")
        .arg(work)
        .current_dir(repo_root())
        .output()
        .unwrap();
    assert!(
        built.status.success(),
        "og live demo failed:\n{}{}",
        String::from_utf8_lossy(&built.stdout),
        String::from_utf8_lossy(&built.stderr)
    );
}

fn call(b: &mut Browser, method: &str, params: &str) -> Value {
    b.call(method, loads(params).unwrap()).unwrap()
}

fn dark(b: &mut Browser) {
    call(
        b,
        "Emulation.setEmulatedMedia",
        r#"{"features":[{"name":"prefers-color-scheme","value":"dark"}]}"#,
    );
}

fn viewport(b: &mut Browser, width: u32, height: u32) {
    call(
        b,
        "Emulation.setDeviceMetricsOverride",
        &format!(r#"{{"width":{width},"height":{height},"deviceScaleFactor":1,"mobile":false}}"#),
    );
}

fn eval(b: &mut Browser, script: &str) -> Value {
    b.eval(script).unwrap_or_else(|e| panic!("{script}: {e}"))
}

fn require(b: &mut Browser, expr: &str) -> Value {
    b.wait(&format!("Boolean({expr})"), WAIT).unwrap_or_else(|| {
        panic!("browser condition never became true: {expr}\n{:?}", b.eval(
            "[location.href, document.documentElement.dataset, document.body.textContent.slice(0,600)]"
        ))
    })
}

fn go(b: &mut Browser, url: &str) {
    b.go("about:blank").unwrap();
    require(b, "location.href === 'about:blank'");
    b.go(url).unwrap();
}

fn run_at(base: &str, sample: &str, suffix: &str) -> String {
    format!("{base}runs/{sample}/{suffix}")
}

fn replay(b: &mut Browser, value: &str) {
    eval(
        b,
        &format!(
            "var s = document.getElementById('replay-slider'); s.value = {value}; s.dispatchEvent(new Event('input')); true"
        ),
    );
    require(
        b,
        "document.documentElement.dataset.pgpFrame === document.getElementById('replay-slider').value",
    );
}

fn css(b: &mut Browser, selector: &str, count: usize) {
    require(
        b,
        &format!(
            "document.querySelectorAll({}).length === {count}",
            js_str(selector)
        ),
    );
}

fn eval_text(b: &mut Browser, expr: &str) -> String {
    if let Value::Str(s) = eval(b, expr) {
        s
    } else {
        panic!("{expr} did not return text");
    }
}

fn assert_no_policy_block(b: &mut Browser) {
    assert_eq!(
        eval(b, "document.documentElement.dataset.pgpBlocked || ''"),
        Value::Str(String::new())
    );
}

fn open_child(b: &mut Browser) {
    require(b, "document.documentElement.dataset.pgpGated === 'true'");
    for _ in 0..10 {
        if eval(b, "!!document.querySelector('[data-opening^=closed]')") != Value::Bool(true) {
            require(b, "!document.documentElement.dataset.pgpGated");
            return;
        }
        // The last click can navigate before CDP replies.
        let _ = b.eval("document.querySelector('[data-opening^=closed]').click(); true");
        if b.wait("!document.documentElement.dataset.pgpGated && document.documentElement.dataset.pgpResults",
            Duration::from_millis(400)).is_some() {
            return;
        }
    }
    panic!("the child report did not open");
}

fn map_level(b: &mut Browser, level: &str) {
    require(
        b,
        &format!(
            "{MAP}?.dataset.pgpMap === 'drawn' && {MAP}.dataset.level === {}",
            js_str(level)
        ),
    );
}

fn replay_and_search(b: &mut Browser, base: &str) {
    go(b, base);
    require(
        b,
        &format!("{READY} && document.documentElement.dataset.pgpPeople === '3'"),
    );
    assert_eq!(
        eval_text(
            b,
            "Array.from(document.querySelectorAll('.person-card a.person-link'), a => a.getAttribute('href')).join(' ')"
        ),
        "runs/DEMO-001/ runs/DEMO-002/ runs/DEMO-003/"
    );
    css(b, ".family-head a[href='runs/DEMO-FAMILY/']", 2);
    go(b, &run_at(base, MAIN, ""));
    require(b, READY);
    replay(b, "0");
    require(
        b,
        "document.querySelectorAll('#map g.box').length > 30 && document.getElementById('run-state').textContent === 'running'",
    );
    assert_eq!(
        eval(b, "document.documentElement.dataset.pgpCharts || ''"),
        Value::Str(String::new())
    );
    let end = "s.max";
    replay(b, "Math.floor(s.max / 2)");
    require(
        b,
        "document.getElementById('content').lastElementChild.id === 'progress' && !document.getElementById('progress-fold').open",
    );
    require(
        b,
        "!document.getElementById('runhero').hidden && !!(document.getElementById('runhero').compareDocumentPosition(document.getElementById('lead-slot')) & Node.DOCUMENT_POSITION_FOLLOWING)",
    );
    require(
        b,
        "document.getElementById('run-state').textContent === 'running'",
    );
    replay(b, end);
    require(b, RESULTS);
    require(b, "document.getElementById('runhero').hidden");
    require(
        b,
        "document.querySelectorAll('#found > li[role=button]').length <= 8",
    );
    require(
        b,
        "document.querySelectorAll('#card-hla .card-text .lines-box > ul.lines > li').length <= 6 && document.querySelector('#card-hla details.lines-more') && document.querySelector('#card-hla .card-margin .caveat .caveat-lead')",
    );
    require(
        b,
        "Number(document.documentElement.dataset.pgpCharts) > 0 && document.getElementById('run-state').textContent === 'finished'",
    );
    require(
        b,
        "document.querySelector('section.plugin[data-plugin=ops]')?.textContent.includes('retried')",
    );
    require(
        b,
        "document.querySelector('section.plugin[data-plugin=recheck]')?.textContent.includes('newer ClinVar')",
    );
    css(
        b,
        "section.plugin[data-plugin=family] a[href='../../runs/DEMO-FAMILY/']",
        1,
    );
    require(
        b,
        "document.querySelector('input[type=search]') && !document.querySelector('input[type=search]').disabled",
    );
    eval(
        b,
        "var i = document.querySelector('input[type=search]'); i.value = 'DEMOC'; i.dispatchEvent(new Event('input')); true",
    );
    require(b, "Number(document.documentElement.dataset.pgpSearch) > 0");
    assert_no_policy_block(b);
    require(
        b,
        "document.querySelector('meta[http-equiv=Content-Security-Policy]').content.includes(\"script-src 'self'\")",
    );

    css(b, "main section.sec[data-journal]", 10);
    css(b, "#rail [data-toc]", 10);
    css(b, "#lead-slot .ledger a[data-level]", 5);
    require(
        b,
        "document.querySelector('#lead-slot .lede')?.textContent.length > 40",
    );
    eval(b, "location.hash = '#s/clinvar'; true");
    require(
        b,
        "document.activeElement.id === 'card-clinvar' && document.querySelector('#card-clinvar .caveat .caveat-lead') && document.querySelector('#card-clinvar figure[data-chart]')",
    );
    let finding = eval_text(
        b,
        "document.querySelector('[data-finding] h3 a[href^=\"finding/\"]')?.getAttribute('href')",
    );
    go(b, &run_at(base, MAIN, &finding));
    require(
        b,
        "document.querySelector('h1.h-display') && document.querySelectorAll('#plan li').length > 0 && document.getElementById('evidence')",
    );
    assert_no_policy_block(b);
}

fn family_view(b: &mut Browser, base: &str) {
    go(b, &run_at(base, FAMILY, ""));
    require(
        b,
        &format!("{READY} && document.querySelector('img.pedigree')?.naturalWidth > 0"),
    );
    css(b, "svg.pedigree-svg a[href='../../runs/DEMO-001/']", 1);
    css(b, "#map, #report-link", 0);
    require(
        b,
        "document.body.textContent.includes('DEMO-001') && document.body.textContent.includes('DEMO-002') && document.body.textContent.includes('DEMO-003')",
    );
    // A parent can see their own report. The child can open their reports
    // independently of the member page that lists them.
    go(b, &run_at(base, PARTNER, ""));
    require(b, RESULTS);
    assert_eq!(
        eval(b, "document.documentElement.dataset.pgpGated || ''"),
        Value::Str(String::new())
    );
    go(b, &run_at(base, CHILD, ""));
    require(b, "document.documentElement.dataset.pgpGated === 'true'");
    css(
        b,
        "#results article[data-finding], #genome-map, .gm-inview, .gm-strip, canvas",
        0,
    );
    open_child(b);
    require(b, RESULTS);
}

fn family_depth(b: &mut Browser, base: &str, site: &Path) {
    fn open(b: &mut Browser, base: &str) {
        go(b, &run_at(base, FAMILY, ""));
        require(b, READY);
        for _ in 0..10 {
            let left = eval(
                b,
                "document.querySelectorAll('[data-opening^=closed]').length",
            );
            let Value::Int(left) = left else {
                panic!("the family gate has no opening count")
            };
            if left == 0 {
                break;
            }
            let _ = b.eval("document.querySelector('[data-opening^=closed]').click(); true");
            require(
                b,
                &format!("document.querySelectorAll('[data-opening^=closed]').length < {left}"),
            );
        }
        require(b, "!document.querySelector('.report-gate')");
    }
    let source = repo_root().join("tests/fixtures/live/demo/DEMO-FAMILY/DEMO-003_segments.tsv");
    let rows = fs::read_to_string(source)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count()
        - 1;
    let paint = "[data-painting='DEMO-003'] .paint-seg";
    open(b, base);
    css(b, "table.rel-grid td:not(.self):not(.none)", 6);
    css(b, paint, rows);
    css(b, "details[data-child='DEMO-003'] table.rows", 3);
    css(b, "section.side-by-side table.rows thead th", 5);

    let path = site.join(format!("runs/{FAMILY}/family.json"));
    let original = fs::read_to_string(&path).unwrap();
    let Value::Object(mut altered) = loads(&original).unwrap() else {
        panic!("family.json is not an object")
    };
    altered.insert("segments".into(), Value::Object(Default::default()));
    fs::write(
        &path,
        pycompat::json::dumps(&Value::Object(altered), &Default::default()).unwrap(),
    )
    .unwrap();
    open(b, base);
    css(b, paint, 0);
    fs::write(path, original).unwrap();
}

fn report_links(b: &mut Browser, base: &str) {
    for run in [MAIN, PARTNER, CHILD] {
        go(b, &run_at(base, run, "report.html"));
        require(
            b,
            "document.querySelector('body > p[role=note]')?.textContent.startsWith('Invented data:')",
        );
        css(b, "a[href='cpsr/']", 1);
        eval(
            b,
            "document.querySelector(\"a[href='cpsr/']\").click(); true",
        );
        require(
            b,
            &format!("document.title === {}", js_str(&format!("CPSR {run}"))),
        );
        require(b, "location.pathname.endsWith('/cpsr/')");
    }
}

fn widget_checks(b: &mut Browser, base: &str, widgets: bool) {
    for (name, checks) in [
        (
            PARTNER,
            &[
                ("#card-expansions [data-widget=str] .w-row", 8),
                ("#card-lpa [data-widget=ci] .w-tick", 1),
                ("#card-chip [data-widget=vaf] .w-tick.w-call", 1),
                ("#card-paralogs [data-widget=copies] li", 7),
                ("#card-clinvar figure[data-chart=genome] rect.mca", 1),
            ][..],
        ),
        (
            CHILD,
            &[
                ("#card-y-chromosome [data-widget=ypath] li", 5),
                ("#card-clinvar figure[data-chart=genome] rect.sv", 2),
            ][..],
        ),
    ] {
        go(b, &run_at(base, name, ""));
        if name == CHILD
            && b.wait(
                "document.documentElement.dataset.pgpGated === 'true'",
                Duration::from_millis(600),
            )
            .is_some()
        {
            open_child(b);
        }
        require(b, RESULTS);
        for (selector, count) in checks {
            let want = if widgets || selector.contains("rect.") {
                *count
            } else {
                0
            };
            css(b, selector, want);
        }
    }
}

fn phone_fits(b: &mut Browser) -> bool {
    eval(
        b,
        "Math.max(document.documentElement.scrollWidth, document.body.getBoundingClientRect().right) <= innerWidth + 1",
    ) == Value::Bool(true)
}

fn phone(b: &mut Browser, base: &str) {
    viewport(b, 390, 844);
    for (url, condition) in [
        (base.to_owned(), READY),
        (run_at(base, MAIN, ""), RESULTS),
        (run_at(base, FAMILY, ""), READY),
        (run_at(base, PARTNER, ""), RESULTS),
        (run_at(base, PARTNER, "finding/str-fmr1/"), READY),
    ] {
        go(b, &url);
        require(b, condition);
        assert!(phone_fits(b), "the phone page overflows at {url}");
    }
    eval(
        b,
        "document.body.style.setProperty('min-width', '900px', 'important'); true",
    );
    assert!(!phone_fits(b), "the wide block control should fail");
    go(b, &run_at(base, PARTNER, ""));
    require(b, RESULTS);
    require(
        b,
        "getComputedStyle(document.getElementById('tabbar')).display === 'grid'",
    );
    require(
        b,
        "document.documentElement.scrollHeight < 8000 && !!document.querySelector('#results details.pfold')",
    );
    eval(b, "location.hash = '#s/roh'; true");
    require(
        b,
        "document.activeElement.id === 'card-roh' && document.getElementById('card-roh').checkVisibility()",
    );
    viewport(b, 1280, 900);
    require(b, "!document.querySelector('#results details.pfold')");
    viewport(b, 390, 844);
    require(
        b,
        "document.querySelector('details[data-card=\"fold-roh\"][open]')",
    );
    call(b, "Emulation.clearDeviceMetricsOverride", "{}");
}

fn map_checks(b: &mut Browser, base: &str) {
    go(b, &run_at(base, MAIN, ""));
    map_level(b, "genome");
    require(
        b,
        "document.getElementById('gm-title').textContent.includes('all 24 chromosomes')",
    );
    let views = [
        ("#chr16", "chromosome"),
        ("#chr11:4.5-6.5M", "region"),
        ("#chr11:5.150-5.350M", "gene"),
        ("#chr11:5225001-5235000", "reads"),
    ];
    for (place, level) in views {
        eval(b, &format!("location.hash = {}; true", js_str(place)));
        map_level(b, level);
        if level != "chromosome" {
            require(b, &format!("Number({MAP}.dataset.genes) > 0"));
        }
    }
    eval(
        b,
        "document.querySelector('.gm-canvas').dispatchEvent(new KeyboardEvent('keydown', {key:'Escape', bubbles:true})); true",
    );
    require(
        b,
        "document.getElementById('genome-map').dataset.level !== 'reads'",
    );
    go(b, &run_at(base, MAIN, "#chr16:2.075-2.100M"));
    require(
        b,
        &format!("{MAP}?.dataset.pgpMap === 'drawn' && Number({MAP}.dataset.genes) > 0"),
    );
    require(b, &format!("{MAP}.textContent.includes('DEMOC')"));
    let source = call(
        b,
        "Page.addScriptToEvaluateOnNewDocument",
        r#"{"source":"delete window.DecompressionStream"}"#,
    );
    go(b, &run_at(base, MAIN, ""));
    require(
        b,
        &format!(
            "{MAP}?.dataset.pgpMap === 'nodecompress' && document.querySelectorAll('.gm-fallback table tbody tr').length > 0"
        ),
    );
    if let Value::Object(o) = source
        && let Some(Value::Str(id)) = o.get("identifier")
    {
        call(
            b,
            "Page.removeScriptToEvaluateOnNewDocument",
            &format!(r#"{{"identifier":{}}}"#, js_str(id)),
        );
    } else {
        panic!("Chrome did not register the missing-decompressor control");
    }
}

fn karyogram_marks(b: &mut Browser, base: &str) {
    let fixture =
        fs::read_to_string(repo_root().join("tests/fixtures/map/karyogram-marks.tsv")).unwrap();
    for person in [MAIN, PARTNER] {
        go(b, &run_at(base, person, ""));
        map_level(b, "genome");
        let mut count = 0;
        for row in fixture.lines().filter(|row| row.starts_with(person)) {
            let parts: Vec<_> = row.split('\t').collect();
            let [_, id, chrom, mb, level] = parts.as_slice() else {
                panic!("bad karyogram fixture row")
            };
            let bp = mb.parse::<f64>().unwrap() * 1e6;
            require(
                b,
                &format!(
                    "[...document.querySelectorAll('.gm-marks a')].some(a => a.dataset.id === {} && a.dataset.chrom === {} && a.dataset.level === {} && Math.abs(Number(a.dataset.s) - {bp}) <= 10000)",
                    js_str(id),
                    js_str(&format!("chr{chrom}")),
                    js_str(level)
                ),
            );
            count += 1;
        }
        assert!(count > 0, "{person} has no karyogram fixture rows");
    }
}

fn map_search(b: &mut Browser, base: &str) {
    fn search(b: &mut Browser, query: &str, condition: &str) {
        eval(
            b,
            &format!(
                "var i = document.getElementById('gm-find'); i.focus(); i.value = {}; i.dispatchEvent(new Event('input')); true",
                js_str(query)
            ),
        );
        require(b, condition);
    }
    go(b, &run_at(base, MAIN, ""));
    map_level(b, "genome");
    eval(
        b,
        "document.body.dispatchEvent(new KeyboardEvent('keydown', {key:'/', bubbles:true})); true",
    );
    require(b, "document.activeElement.id === 'gm-find'");
    for (query, kind, chrom) in [
        ("HBB", "gene", "chr11"),
        ("11p15.4", "band", "chr11"),
        ("DEMOC", "finding", "chr16"),
    ] {
        let condition = format!(
            "[...document.querySelectorAll('#gm-find-list li[data-kind={kind}]')].some(el => el.textContent.startsWith({}))",
            js_str(query)
        );
        search(b, query, &condition);
        eval(
            b,
            &format!(
                "[...document.querySelectorAll('#gm-find-list li[data-kind={kind}]')].find(el => el.textContent.startsWith({})).dispatchEvent(new MouseEvent('mousedown', {{bubbles:true, cancelable:true}})); true",
                js_str(query)
            ),
        );
        require(b, &format!("{MAP}.dataset.chrom === {}", js_str(chrom)));
        if kind == "finding" {
            require(
                b,
                "document.querySelector('.gm-card').textContent.includes('DEMOC')",
            );
        }
    }
    search(
        b,
        "APOE",
        "document.querySelector('#gm-find-list li[data-kind=gene]')?.textContent.startsWith('APOE')",
    );
    css(b, "#gm-find-list li[data-kind=finding]", 0);
}

const IN_VIEW_MATCHES_CANVAS: &str = r#"(() => {
    const map = document.getElementById('genome-map');
    const table = map.querySelector('.gm-inview');
    const canvas = map.querySelector('.gm-canvas');
    const show = table.querySelector('.gm-showall');
    if (show?.textContent === 'Show all') show.click();
    const rows = [...table.querySelectorAll('tbody tr[data-key]')].map(row => row.dataset.key);
    const drawn = canvas.pgpDrawn || [];
    const count = rows.length + (rows.length === 1 ? ' feature' : ' features');
    return new Set(rows).size === rows.length
        && table.querySelector('.gm-inview-n').textContent === count
        && JSON.stringify(rows.sort()) === JSON.stringify(drawn.slice().sort());
})()"#;

fn map_table(b: &mut Browser, base: &str) {
    go(b, &run_at(base, MAIN, ""));
    map_level(b, "genome");
    require(b, IN_VIEW_MATCHES_CANVAS);
    require(b, "document.querySelector('.gm-strip').hidden");
    require(
        b,
        "document.querySelector('.gm-inview-dense').textContent.includes('shorter runs of homozygosity')",
    );
    let marks = eval(
        b,
        "document.querySelector('.gm-canvas').pgpDrawn.filter(key => key.startsWith('findings\\t')).sort()",
    );
    assert!(
        matches!(&marks, Value::Array(v) if !v.is_empty()),
        "genome has no open findings"
    );
    eval(b, "location.hash = '#chr19'; true");
    map_level(b, "chromosome");
    require(b, IN_VIEW_MATCHES_CANVAS);
    require(b, "document.querySelector('.gm-strip').clientHeight === 34");
    assert_eq!(
        eval(
            b,
            "document.querySelector('.gm-strip').pgpDrawn.slice().sort()"
        ),
        marks
    );
    require(
        b,
        r#"![...document.querySelectorAll('.gm-table tr[data-key]')]
        .some(row => row.dataset.key.includes('\tapoe'))
        && !document.querySelector('.gm-strip').pgpDrawn.some(key => key.includes('\tapoe'))"#,
    );
    eval(b, "location.hash = '#chr22:40.5-43.5M'; true");
    map_level(b, "region");
    require(b, &format!("Number({MAP}.dataset.genes) > 0"));
    require(b, IN_VIEW_MATCHES_CANVAS);
    require(
        b,
        r#"(() => {
        const live = document.querySelector('.gm-inview-live');
        const table = document.querySelector('.gm-inview');
        return live && !live.closest('details') && table.open;
    })()"#,
    );
    eval(
        b,
        "document.body.dispatchEvent(new KeyboardEvent('keydown', {key:'t', bubbles:true})); true",
    );
    require(b, "!document.querySelector('.gm-inview').open");
    eval(
        b,
        "document.body.dispatchEvent(new KeyboardEvent('keydown', {key:'t', bubbles:true})); true",
    );
    require(b, "document.querySelector('.gm-inview').open");
    eval(
        b,
        "document.querySelector('.gm-table tbody tr[data-key]').remove(); true",
    );
    assert_eq!(
        eval(b, IN_VIEW_MATCHES_CANVAS),
        Value::Bool(false),
        "missing row control should fail"
    );
}

const PREDICATE_ROWS: &str = r#"(async function () {
    const out = {}, genes = [], {LEVELS, drawsAt} = window.PGPMap;
    async function gz(url) {
        const response = await fetch(url, {cache: 'no-store'});
        if (!response.ok) throw new Error(url + ': ' + response.status);
        return new Response(response.body.pipeThrough(new DecompressionStream('gzip'))).json();
    }
    function row(layer, chrom, start, end, id, data) {
        return [layer, chrom, start, end, id, ...LEVELS.map(level => +drawsAt(layer, data, level))].join('\t');
    }
    function each(doc, fn) {
        let start = 0;
        (doc.ds || []).forEach((delta, i) => { start += delta; fn(start, start + doc.l[i], i); });
    }
    const chroms = Array.from({length:22}, (_, i) => 'chr' + (i + 1)).concat('chrX','chrY');
    const geneFiles = await Promise.all(chroms.map(chrom => gz('../../map/genes/' + chrom + '.json')));
    geneFiles.forEach(({g}, i) => each(g, (s,e,j) => {
        genes.push(row('genes', chroms[i], s, e, g.name[j], {name:g.name[j]}));
    }));
    for (const run of ['DEMO-001','DEMO-002','DEMO-003']) {
        const {per} = await gz('../' + run + '/map.json');
        const {findings:f, roh} = per[run];
        const rows = f.id.map((id, i) => row('findings', f.chrom[i], f.s[i], f.e[i], id,
            {id, level:f.level[i], s:f.s[i], e:f.e[i]}));
        for (const [chrom, blocks] of Object.entries(roh || {})) {
            each(blocks, (s,e,i) => rows.push(row('roh', chrom, s, e, '', {long:blocks.long[i]})));
        }
        out[run] = rows.concat(genes).sort();
    }
    return out;
})()"#;

fn map_command(og: &Path, site: &Path, work: &Path, run: &str, predicate: bool) -> Vec<u8> {
    let results = fs::read_to_string(site.join(format!("runs/{run}/results.json"))).unwrap();
    let mut command = Command::new(og);
    command
        .args(["map", "build", "--root"])
        .arg(work)
        .args(["--run", run, "--assets"])
        .arg(repo_root().join("assets"))
        .arg("--cache-dir")
        .arg(work.join(".owngenome-live/map"));
    if predicate {
        command.arg("--predicate-table");
    }
    let mut child = command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    write!(child.stdin.take().unwrap(), "{{\"results\":{results}}}").unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    result.stdout
}

fn map_predicate(b: &mut Browser, base: &str, og: &Path, site: &Path, work: &Path) {
    let map = site.join(format!("runs/{MAIN}/map.json"));
    let original = fs::read(&map).unwrap();
    let table = work.join(format!("{MAIN}/vcf/{MAIN}_roh.txt"));
    let mut file = fs::OpenOptions::new().append(true).open(&table).unwrap();
    writeln!(
        file,
        "RG\t{MAIN}\tchr2\t60000001\t67000000\t7000000\t5000\t80.0"
    )
    .unwrap();
    drop(file);
    fs::write(&map, map_command(og, site, work, MAIN, false)).unwrap();
    go(b, &run_at(base, MAIN, ""));
    map_level(b, "genome");
    let response = call(
        b,
        "Runtime.evaluate",
        &format!(
            r#"{{"expression":{},"awaitPromise":true,"returnByValue":true}}"#,
            js_str(PREDICATE_ROWS)
        ),
    );
    let Some(Value::Object(javascript)) = prop(&response, "result").and_then(|v| prop(v, "value"))
    else {
        panic!("the map predicate did not return its rows: {response:?}");
    };
    let mut seen = std::collections::BTreeSet::new();
    for run in [MAIN, PARTNER, CHILD] {
        let output = String::from_utf8(map_command(og, site, work, run, true)).unwrap();
        let mut lines = output.lines();
        assert_eq!(
            lines.next(),
            Some("layer\tchrom\ts\te\tid\tgenome\tchromosome\tregion\tgene\treads")
        );
        let mut rust: Vec<_> = lines.map(str::to_owned).collect();
        rust.sort_unstable();
        let Some(Value::Array(js)) = javascript.get(run) else {
            panic!("{run} has no browser predicate rows");
        };
        let js: Vec<_> = js
            .iter()
            .filter_map(|v| {
                if let Value::Str(v) = v {
                    Some(v.as_str())
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(
            rust.len(),
            js.len(),
            "{run} has a different number of map rows"
        );
        for (row, actual) in rust.iter().zip(js) {
            assert_eq!(row, actual, "different map predicate for {run}");
            let row: Vec<_> = row.split('\t').collect();
            seen.insert((row[0].to_owned(), row[5].to_owned()));
        }
    }
    for (layer, visible) in [
        ("findings", "1"),
        ("roh", "1"),
        ("roh", "0"),
        ("genes", "0"),
    ] {
        assert!(
            seen.contains(&(layer.to_owned(), visible.to_owned())),
            "missing predicate input for {layer}/{visible}"
        );
    }
    fs::write(map, original).unwrap();
}

fn prop<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    if let Value::Object(o) = v {
        o.get(key)
    } else {
        None
    }
}

#[test]
fn replay_reports_and_map_open_under_a_subpath() {
    let Some(chrome) = chrome_or_skip("open the public demo") else {
        return;
    };
    let Some(og) = og() else { return };
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "public-demo");
    let demo = scratch.path().join("demo");
    let work = scratch.path().join("run");
    build(&og, &demo, &work);
    let server = StaticServer::start(scratch.path(), &[("Cache-Control", "no-store")]).unwrap();
    let base = server.url("demo/");
    let mut b = Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();
    dark(&mut b);
    viewport(&mut b, 1440, 900);
    call(&mut b, "Page.enable", "{}");
    replay_and_search(&mut b, &base);
    family_view(&mut b, &base);
    family_depth(&mut b, &base, &demo);
    report_links(&mut b, &base);
    widget_checks(&mut b, &base, true);
    phone(&mut b, &base);
    map_checks(&mut b, &base);
    map_table(&mut b, &base);
    map_search(&mut b, &base);
    karyogram_marks(&mut b, &base);
    map_predicate(&mut b, &base, &og, &demo, &work);
    perf::check(&mut b, &base, &demo);

    let person = demo.join("assets/person.js");
    let original = fs::read_to_string(&person).unwrap();
    let widget = "add(text, widget(c.widget));";
    assert!(
        original.contains(widget),
        "widget control did not patch the script"
    );
    fs::write(&person, original.replace(widget, "")).unwrap();
    widget_checks(&mut b, &base, false);
    fs::write(&person, original).unwrap();
}

fn screenshot(b: &mut Browser, file: &Path) {
    let Value::Object(picture) = call(b, "Page.captureScreenshot", r#"{"format":"png"}"#) else {
        panic!("Chrome returned no screenshot object");
    };
    let Some(Value::Str(encoded)) = picture.get("data") else {
        panic!("Chrome returned no screenshot data");
    };
    let mut decoder = Command::new("base64")
        .arg("--decode")
        .stdin(Stdio::piped())
        .stdout(fs::File::create(file).unwrap())
        .spawn()
        .unwrap();
    decoder
        .stdin
        .take()
        .unwrap()
        .write_all(encoded.as_bytes())
        .unwrap();
    assert!(decoder.wait().unwrap().success());
    assert!(fs::read(file).unwrap().starts_with(b"\x89PNG\r\n\x1a\n"));
}

#[test]
#[ignore = "capture the Docs workflow artifact from an already built demo"]
fn capture_map_screens() {
    let Some(chrome) = chrome_or_skip("take the demo map screenshots") else {
        return;
    };
    let site = PathBuf::from(std::env::var_os("OG_DEMO_SITE").expect("set OG_DEMO_SITE"));
    let dest = PathBuf::from(std::env::var_os("OG_DEMO_SCREENS").expect("set OG_DEMO_SCREENS"));
    fs::create_dir_all(&dest).unwrap();
    let scratch = Scratch::new(Path::new(env!("CARGO_TARGET_TMPDIR")), "map-screens");
    let srv = StaticServer::start(&site, &[]).unwrap();
    let mut b = Browser::launch(&chrome, &scratch.path().join("chrome")).unwrap();
    call(&mut b, "Page.enable", "{}");
    dark(&mut b);
    for (width, height) in [(1440, 900), (390, 844)] {
        viewport(&mut b, width, height);
        for (sample, hash, name, level) in [
            (MAIN, "", "genome", "genome"),
            (MAIN, "#chr16", "chromosome", "chromosome"),
            (MAIN, "#chr11:4.5-6.5M", "region", "region"),
            (MAIN, "#chr11:5.150-5.350M", "gene", "gene"),
            (CHILD, "", "closed", "gate"),
            (MAIN, "#chr22:40.5-43.5M", "inview", "region"),
        ] {
            go(&mut b, &srv.url(&format!("runs/{sample}/{hash}")));
            if level == "gate" {
                require(
                    &mut b,
                    "document.documentElement.dataset.pgpGated === 'true'",
                );
                eval(
                    &mut b,
                    "document.getElementById('report-gate').scrollIntoView({block:'start'}); true",
                );
            } else {
                map_level(&mut b, level);
                if matches!(level, "region" | "gene") {
                    require(&mut b, &format!("Number({MAP}.dataset.genes) > 0"));
                }
                let script = if name == "inview" {
                    "var t = document.querySelector('.gm-inview'); t.open = true; t.scrollIntoView({block:'start'}); window.scrollBy(0,-80); true"
                } else {
                    "document.getElementById('genome-map').scrollIntoView({block:'start'}); true"
                };
                eval(&mut b, script);
            }
            let file = dest.join(format!("map-{name}-{width}x{height}.png"));
            screenshot(&mut b, &file);
            println!("{}", file.display());
        }
    }
}
