//! The CSS of the HTML report: [`CSS`] byte for byte what `bin/render_report.py`
//! writes into the page's `<style>` element: a dark page, a light printout,
//! and print rules that keep one section or chapter when the address names it;
//! then [`READING`], the rules for a phone and a keyboard the Python never had.

/// The `<style>` body, from its leading newline to its closing brace and newline.
pub const CSS: &str = r#"
  :root { --bg: #15171c; --card: #1e2128; --fg: #e3e5ea; --muted: #a3a9b4; --line: #2e333c; --head: #262a32;
          --accent: #8cbcf0; --hdr: linear-gradient(135deg, #12324d 0%, #1f5f8f 100%); --hdr-fg: #f2f6fa;
          --green-bg: #1d3b2a; --green-fg: #8ee0ae; --yellow-bg: #3f3618; --yellow-fg: #f2d77c;
          --red-bg: #4a2424; --red-fg: #f5a9a2; --gray-bg: #2d3139; --gray-fg: #b6bbc4;
          --state-bg: #1f3550; --state-fg: #a9cdf5; --stale-bg: #3c2e16; --stale-fg: #f4ca8c;
          --caveat-bg: #232a35; --disc-bg: #3a3214; --disc-line: #a8861c; color-scheme: dark; }
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         background: var(--bg); color: var(--fg); line-height: 1.6; }
  a { color: var(--accent); }
  .layout { display: grid; grid-template-columns: 250px minmax(0, 1fr); gap: 24px; max-width: 1400px;
            margin: 0 auto; padding: 20px; }
  .toc { position: sticky; top: 0; align-self: start; max-height: 100vh; overflow-y: auto; padding: 8px 0;
         font-size: 13px; }
  .toc .toc-title { font-size: 15px; color: var(--accent); margin-bottom: 8px; }
  .toc ol { list-style: none; }
  .toc > ol > li { margin-bottom: 10px; }
  .toc > ol > li > a { font-weight: 600; }
  .toc ol ol { margin: 2px 0 0 4px; }
  .toc a { color: var(--fg); text-decoration: none; }
  .toc a:hover { color: var(--accent); }
  .dot { display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 6px; }
  .dot-ok { background: #5b9ee6; } .dot-stale, .dot-incomplete { background: #e0b94a; }
  .dot-failed, .dot-part { background: #e2685d; } .dot-missing { background: #6b717c; }
  .header { background: var(--hdr); color: var(--hdr-fg); padding: 26px 30px; border-radius: 12px; margin-bottom: 16px; }
  .header h1 { font-size: 28px; margin-bottom: 8px; }
  .header .meta { opacity: 0.9; font-size: 14px; }
  .print-hint { font-size: 13px; color: var(--muted); margin-bottom: 20px; }
  .chapter { margin-bottom: 32px; }
  .chapter-title { font-size: 22px; color: var(--accent); }
  .question { color: var(--muted); margin-bottom: 14px; }
  .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 20px; }
  .sec { display: flex; flex-direction: column; gap: 20px; scroll-margin-top: 12px; }
  .sec:target .card { outline: 2px solid var(--accent); }
  .card { background: var(--card); border-radius: 10px; padding: 24px; box-shadow: 0 2px 8px rgba(0,0,0,0.3); }
  .card h2 { font-size: 18px; color: var(--accent); margin-bottom: 8px;
             padding-bottom: 8px; border-bottom: 2px solid var(--line); }
  .card h3 { font-size: 15px; margin: 12px 0 6px; }
  .sec-meta { font-size: 13px; color: var(--muted); margin-bottom: 8px; }
  .sec-meta .anchor { color: var(--muted); }
  .caveat { font-size: 13px; background: var(--caveat-bg); border-left: 3px solid var(--accent);
            padding: 6px 10px; margin-bottom: 10px; border-radius: 4px; }
  .note { font-size: 13px; color: var(--muted); margin-top: 8px; }
  .files { font-size: 13px; color: var(--muted); }
  .offsite-note { font-size: 12px; color: var(--muted); font-weight: normal; }
  .stat { display: flex; justify-content: space-between; padding: 8px 0;
          border-bottom: 1px solid var(--line); gap: 12px; }
  .stat:last-child { border-bottom: none; }
  .stat .label { color: var(--muted); }
  .stat .value { font-weight: 600; text-align: right; }
  .badge { display: inline-block; padding: 2px 10px; border-radius: 12px;
           font-size: 13px; font-weight: 600; }
  .badge-green { background: var(--green-bg); color: var(--green-fg); }
  .badge-yellow { background: var(--yellow-bg); color: var(--yellow-fg); }
  .badge-red { background: var(--red-bg); color: var(--red-fg); }
  .badge-gray { background: var(--gray-bg); color: var(--gray-fg); }
  .badge-state { background: var(--state-bg); color: var(--state-fg); }
  .stale { background: var(--stale-bg); color: var(--stale-fg); border-radius: 6px; padding: 6px 10px;
           font-size: 13px; margin-bottom: 10px; }
  table { width: 100%; border-collapse: collapse; font-size: 14px; }
  th, td { padding: 8px 12px; text-align: left; border-bottom: 1px solid var(--line); }
  th { background: var(--head); font-weight: 600; color: var(--muted); }
  .full-width { grid-column: 1 / -1; }
  figure.chart { margin: 4px 0 12px; } figure.chart svg { display: block; max-width: 760px; }
  figure.chart a:focus, figure.chart a:hover { opacity: 0.8; }
  ul.plain { margin-left: 18px; }
  .disclaimer { background: var(--disc-bg); border: 1px solid var(--disc-line); border-radius: 8px;
                padding: 16px; margin-top: 24px; font-size: 14px; }
  .footer { text-align: center; color: var(--muted); font-size: 13px; margin-top: 24px; padding: 16px; }
  .footer table { font-size: 12px; color: var(--muted); margin-top: 8px; }
  @media (max-width: 900px) { .layout { display: block; } .toc { position: static; max-height: none; } }
  @media (max-width: 700px) { .grid { grid-template-columns: 1fr; } }
  /* Print is light. With a section or a chapter in the address (#id, from
     the contents), only it prints, with the header, the list of what is not
     assessed (.always) and the disclaimer. A browser without :has() prints
     everything. */
  @media print {
    :root { --bg: #fff; --card: #fff; --fg: #111; --muted: #444; --line: #ccc; --head: #f0f0f0; --accent: #1a5276;
            --hdr: #fff; --hdr-fg: #111; --green-bg: #d5f5e3; --green-fg: #196f3d; --yellow-bg: #fef9e7;
            --yellow-fg: #7d6608; --red-bg: #fadbd8; --red-fg: #922b21; --gray-bg: #eee; --gray-fg: #555;
            --state-bg: #e3eefa; --state-fg: #1a4a7a; --stale-bg: #fdebd0; --stale-fg: #7e5109;
            --caveat-bg: #f4f6f8; --disc-bg: #fff3cd; --disc-line: #c9a227; color-scheme: light; }
    .toc, .print-hint, .sec-meta .anchor { display: none; }
    .layout { display: block; padding: 0; max-width: none; }
    .header { border: 1px solid #999; padding: 12px 16px; }
    .card { box-shadow: none; border: 1px solid #ccc; break-inside: avoid; }
    .sec:target .card { outline: none; }
    .grid { display: block; } .sec { margin-bottom: 16px; }
    body:has(main :target) .chapter:not(:target):not(:has(:target)):not(:has(.always)) { display: none; }
    body:has(.sec:target) .sec:not(:target):not(.always) { display: none; }
    body:has(.chapter:target) .chapter:not(:target) .sec:not(.always) { display: none; }
    body:has(main :target) .chapter:not(:target):not(:has(:target)) > .chapter-title,
    body:has(main :target) .chapter:not(:target):not(:has(:target)) > .question { display: none; }
  }
"#;

/// Rules after [`CSS`] in the same `<style>` element, for reading on a phone
/// and from the keyboard: the skip link shows only with the focus; the
/// research-only line under the title; sections, cards and table cells shrink
/// to the screen and wrap long words (paths, coordinates); each table
/// scrolls inside its region; the contents fold under 900 px and show
/// without their summary line above it where the browser can style
/// `::details-content` (elsewhere the summary opens them); a printout keeps
/// whole tables and leaves out the skip link and the contents.
pub const READING: &str = r#"  .skip { position: absolute; left: 8px; top: -60px; z-index: 10; padding: 8px 14px; border-radius: 6px;
           background: var(--card); color: var(--accent); font-weight: 600; }
  .skip:focus { top: 8px; }
  main:focus { outline: none; }
  .header .research-only { font-weight: 600; margin-bottom: 8px; }
  .sec, .card { min-width: 0; overflow-wrap: anywhere; }
  td { overflow-wrap: anywhere; }
  .tscroll { overflow-x: auto; max-width: 100%; }
  .tscroll:focus-visible { outline: 2px solid var(--accent); }
  .toc-fold > summary { cursor: pointer; color: var(--accent); font-weight: 600; padding: 8px 0; }
  @media (max-width: 900px) { .toc-fold .toc-title { display: none; } }
  @media (min-width: 901px) { .toc-fold { position: sticky; top: 0; align-self: start; max-height: 100vh; overflow-y: auto; } }
  @supports selector(::details-content) {
    @media (min-width: 901px) {
      .toc-fold > summary { display: none; }
      .toc-fold::details-content { content-visibility: visible; display: block; }
    }
  }
  @media print { .toc-fold, .skip { display: none; } .tscroll { overflow: visible; } }
"#;
