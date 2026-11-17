// The family page keeps a child's closed de novo and compound-het rows off the page (assets/live/core.js
// projectFamily, assets/live/family-page.js). Node only: core.js and family-page.js run in a vm context that
// stands in for the browser window, with a small DOM, on an invented family (DEMO-FAMILY names, no real data).
//
// The family: DEMO-001 and DEMO-002, their son DEMO-003 (one result to act on, in KCNH2), his partner DEMO-004
// and their daughter DEMO-005 (no result to act on). Checks:
//   1. While DEMO-003's result is closed, neither his trio rows nor DEMO-005's (read against his genome) reach
//      the family view or the page, and the page is the same with or without those rows: no count, no heading.
//   2. Opened on the page with its button, both children's rows show.
//   3. With no result to act on, every row shows.
// Controls: the same checks on a core.js with the trio rule cut out must fail.
//
// Run: node tests/live/family_hide_closed.js [DIR], DIR holding assets/live (default: this repo).
'use strict';
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const REPO = path.resolve(process.argv[2] || path.join(__dirname, '..', '..'));
const CORE = fs.readFileSync(path.join(REPO, 'assets/live/core.js'), 'utf8');
const PAGE = fs.readFileSync(path.join(REPO, 'assets/live/family-page.js'), 'utf8');
let fails = 0;
function check(desc, ok, detail) {
  console.log((ok ? 'PASS ' : 'FAIL ') + desc + (ok || !detail ? '' : ' -- ' + String(detail).slice(0, 600)));
  if (!ok) { fails++; }
}

// --- a small DOM: enough for el/add/svg, textContent, attributes, details and buttons ------------------------
class Node {
  constructor(tag) {
    this.tag = tag; this.attrs = {}; this.kids = []; this.text = ''; this.listeners = {}; this.parent = null;
    this.style = { setProperty() {} };
    this.classList = { add: (c) => { this.attrs.class = ((this.attrs.class || '') + ' ' + c).trim(); }, contains: () => false };
  }
  setAttribute(k, v) { this.attrs[k] = String(v); }
  getAttribute(k) { return k in this.attrs ? this.attrs[k] : null; }
  removeAttribute(k) { delete this.attrs[k]; }
  get attributes() { return Object.keys(this.attrs).map((name) => ({ name, value: this.attrs[name] })); }
  appendChild(c) { if (c.parent) { c.remove(); } c.parent = this; this.kids.push(c); return c; }
  remove() { if (this.parent) { this.parent.kids = this.parent.kids.filter((k) => k !== this); this.parent = null; } }
  replaceWith(o) { const p = this.parent, i = p.kids.indexOf(this); o.parent = p; p.kids[i] = o; this.parent = null; }
  addEventListener(k, f) { (this.listeners[k] = this.listeners[k] || []).push(f); }
  click() { (this.listeners.click || []).forEach((f) => f({ preventDefault() {} })); }
  get textContent() { return this.tag === '#text' ? this.text : this.kids.map((k) => k.textContent).join(''); }
  set textContent(t) { this.kids = []; if (this.tag === '#text') { this.text = String(t); } else if (t !== '') { const n = new Node('#text'); n.text = String(t); this.appendChild(n); } }
  all() { return [this].concat(...this.kids.map((k) => k.all())); }
  // tag[attr] or [attr] or [attr=value]: what the views ask for.
  querySelectorAll(sel) {
    const m = /^([a-z0-9]*)(?:\[([a-z-]+)(?:=([^\]]+))?\])?$/.exec(sel);
    if (!m) { return []; }
    return this.all().slice(1).filter((n) => (!m[1] || n.tag === m[1]) && (!m[2] || (m[2] in n.attrs && (m[3] === undefined || n.attrs[m[2]] === m[3]))));
  }
  querySelector(sel) { return this.querySelectorAll(sel)[0] || null; }
}

// --- one page load: the window, the server's answers, core.js then family-page.js --------------------------
function load(core, files) {
  const html = new Node('html'), body = new Node('body'), main = new Node('main');
  main.attrs.id = 'family';
  body.classList = { contains: (c) => c === 'page-family', add() {} };
  body.appendChild(main); html.appendChild(body);
  const location = { href: 'http://demo.invalid/runs/DEMO-FAMILY/', pathname: '/runs/DEMO-FAMILY/', hash: '', search: '', replace() {} };
  const store = {};
  const win = {
    document: { documentElement: html, body, readyState: 'complete', hidden: false, title: '', addEventListener() {},
      createElement: (t) => new Node(t), createElementNS: (ns, t) => new Node(t),
      createTextNode: (t) => { const n = new Node('#text'); n.text = String(t); return n; },
      getElementById: (id) => html.all().filter((n) => n.attrs.id === id)[0] || null },
    location, name: '', scrollY: 0, scrollTo() {}, setInterval() { return 0; }, addEventListener() {},
    history: { replaceState() {} },
    sessionStorage: { getItem: (k) => store[k] || null, setItem: (k, v) => { store[k] = String(v); }, removeItem: (k) => { delete store[k]; } },
    EventSource: function () {}, URL, Response, TextEncoder, crypto: globalThis.crypto, setTimeout, clearTimeout, console,
    fetch: (u) => {
      const p = new URL(typeof u === 'string' ? u : u.url, location.href).pathname;
      const body = files[p];
      return Promise.resolve(body === undefined ? new Response('not found', { status: 404 }) : new Response(JSON.stringify(body), { status: 200 }));
    }
  };
  win.window = win; win.self = win;
  vm.createContext(win);
  vm.runInContext(core, win, { filename: 'core.js' });
  vm.runInContext(PAGE, win, { filename: 'family-page.js' });
  return win;
}
function tick() { return new Promise((r) => setTimeout(r, 5)); }
async function settle(win, want) {
  for (let i = 0; i < 400; i++) {
    if (win.document.documentElement.attrs['data-pgp-status'] === 'ready' && (!want || want())) { return true; }
    await tick();
  }
  return false;
}
function pageOf(win) {
  const main = win.document.getElementById('family');
  return { text: main.textContent, headings: main.querySelectorAll('h3').map((h) => h.textContent), main };
}
// What the family view itself receives (core.js projects family.json for every reader of it).
async function viewData(win) {
  const r = await win.fetch('/runs/DEMO-FAMILY/family.json');
  return r.json();
}

// --- the invented family ---------------------------------------------------------------------------------
const DN = (child, pos) => ({ child, chrom: 'chr3', pos: String(pos), ref: 'C', alt: 'T', child_gt: '0/1', father_gt: '0/0', mother_gt: '0/0',
  child_dp: '34', child_gq: '52', child_ab: '0.47', quality: 'high' });
const CH = (child, gene, a, b) => ({ child, gene, variant_a: a, variant_b: b, a_from: 'father', b_from: 'mother', phase: 'in trans' });
function family(rows) {
  const members = [
    { sample: 'DEMO-001', sex: 'male', father: '', mother: '', partner: 'DEMO-002', sequenced: true },
    { sample: 'DEMO-002', sex: 'female', father: '', mother: '', partner: 'DEMO-001', sequenced: true },
    { sample: 'DEMO-003', sex: 'male', father: 'DEMO-001', mother: 'DEMO-002', partner: 'DEMO-004', sequenced: true },
    { sample: 'DEMO-004', sex: 'female', father: '', mother: '', partner: 'DEMO-003', sequenced: true },
    { sample: 'DEMO-005', sex: 'female', father: 'DEMO-003', mother: 'DEMO-004', partner: '', sequenced: true }];
  const assessed = { de_novo: { assessed: true, why: 'both parents are sequenced' }, compound_het: { assessed: true, why: 'a parent is sequenced' } };
  const run = (s, child) => ({ url: '/runs/' + s + '/', state: 'finished', planned: 10, done: 10, attention: child ? null : 1, results: true, child });
  return {
    state: { kind: 'family', name: 'DEMO-FAMILY', pedigree: null, members,
      children: [{ child: 'DEMO-003', father: 'DEMO-001', mother: 'DEMO-002' }, { child: 'DEMO-005', father: 'DEMO-003', mother: 'DEMO-004' }],
      couples: [], plan: [], reports: [] },
    layout: null, relationships: [], ibd_pairs: [], couples: {}, segments: {},
    de_novo: rows.de_novo, compound_het: rows.compound_het,
    trio_files: { 'de_novo.tsv': true, 'compound_het.tsv': true }, trio_file_state: { 'de_novo.tsv': 'here', 'compound_het.tsv': 'here' },
    trio_status: { 'DEMO-003': assessed, 'DEMO-005': assessed }, side_by_side: null,
    runs: { 'DEMO-001': run('DEMO-001', false), 'DEMO-002': run('DEMO-002', false), 'DEMO-003': run('DEMO-003', true),
      'DEMO-004': run('DEMO-004', false), 'DEMO-005': run('DEMO-005', true) } };
}
const ROWS = {
  de_novo: [DN('DEMO-003', 48211904), DN('DEMO-005', 51122233)],
  compound_het: [CH('DEMO-003', 'KCNH2', 'chr7:150960000:C>T', 'chr7:150962400:G>A'), CH('DEMO-005', 'G0187', 'chr7:94371020:G>A', 'chr7:94398114:T>C')]
};
// The values a row of each child would put on the page.
const SEEN = { 'DEMO-003': ['48211904', 'KCNH2', '150960000'], 'DEMO-005': ['51122233', 'G0187', '94371020'] };
const ACT = { id: 'act-kcnh2', level: 'act', gene: 'KCNH2', title: 'A change in KCNH2', detail: 'chr7:150960000 C>T', text: 'Invented.', plan: ['Invented.'] };
function files(rows, act) {
  const out = { '/runs/DEMO-FAMILY/family.json': family(rows) };
  out['/runs/DEMO-003/results.json'] = { findings: (act ? [ACT] : []).concat([{ id: 'look-1', level: 'review', gene: 'G0001', title: 'Invented' }]) };
  out['/runs/DEMO-005/results.json'] = { findings: [{ id: 'look-2', level: 'review', gene: 'G0002', title: 'Invented' }] };
  return out;
}
const NONE = { de_novo: [], compound_het: [] };
function rowsOf(d, child) { return (d.de_novo || []).concat(d.compound_het || []).filter((r) => r.child === child).length; }
function shows(text, child) { return SEEN[child].filter((v) => text.indexOf(v) >= 0); }

// --- the checks, on one core.js ------------------------------------------------------------------------------
async function closedChecks(core, label) {
  const ok = [];
  const win = load(core, files(ROWS, true));
  if (!(await settle(win))) { return [false]; }
  const p = pageOf(win), d = await viewData(win);
  ok.push(['closed: no row of DEMO-003 or DEMO-005 reaches the family view (' + label + ')',
    rowsOf(d, 'DEMO-003') === 0 && rowsOf(d, 'DEMO-005') === 0, JSON.stringify([d.de_novo, d.compound_het])]);
  ok.push(['closed: no value of their rows is on the page (' + label + ')',
    !shows(p.text, 'DEMO-003').length && !shows(p.text, 'DEMO-005').length, shows(p.text, 'DEMO-003').concat(shows(p.text, 'DEMO-005'))]);
  ok.push(['closed: no trio heading, count or row button (' + label + ')',
    !p.headings.some((h) => /De novo|Compound-het/.test(h)) && !/Show all|more rows/.test(p.text), p.headings]);
  // The page with the rows and the page without any row are the same: nothing on it says a row is there.
  const empty = load(core, files(NONE, true));
  await settle(empty);
  ok.push(['closed: the page is the same with or without their rows (' + label + ')', pageOf(empty).text === p.text,
    pageOf(empty).text.length + ' vs ' + p.text.length]);
  ok.push(['closed: DEMO-005 keeps her entry and one fixed closed line (' + label + ')',
    p.main.querySelectorAll('[data-trio=closed]').length === 1 && /DEMO-003’s report is closed/.test(p.text), p.text.slice(-400)]);
  return ok;
}
async function openedChecks(core) {
  const ok = [];
  const win = load(core, files(ROWS, true));
  await settle(win);
  const buttons = pageOf(win).main.querySelectorAll('[data-opening]');
  ok.push(['opened: DEMO-003 has one button', buttons.length === 1, buttons.length]);
  buttons.forEach((b) => b.click());
  const done = await settle(win, () => shows(pageOf(win).text, 'DEMO-003').length === 3);
  const p = pageOf(win);
  ok.push(['opened: both children\'s rows show once his result is opened', done && shows(p.text, 'DEMO-005').length === 3,
    shows(p.text, 'DEMO-003').concat(shows(p.text, 'DEMO-005'))]);
  ok.push(['opened: the closed line is gone', !p.main.querySelectorAll('[data-trio=closed]').length]);
  return ok;
}
async function noActChecks(core) {
  const win = load(core, files(ROWS, false));
  await settle(win);
  const p = pageOf(win);
  return [['no result to act on: every row shows', shows(p.text, 'DEMO-003').length === 3 && shows(p.text, 'DEMO-005').length === 3 &&
    p.headings.filter((h) => /De novo|Compound-het/.test(h)).length === 4, p.headings]];
}

(async () => {
  for (const c of await closedChecks(CORE, 'core.js')) { check(c[0], c[1], c[2]); }
  for (const c of await openedChecks(CORE)) { check(c[0], c[1], c[2]); }
  for (const c of await noActChecks(CORE)) { check(c[0], c[1], c[2]); }
  // Controls: the trio rule cut out of core.js (rows kept whoever is gated) must fail the closed checks.
  const RULE = 'function trioShut(s) { return gated(s) || parentsOf(s).some(function (p) { return p && gated(p); }); }';
  check('control: core.js holds the trio rule the control cuts out', CORE.indexOf(RULE) >= 0);
  if (CORE.indexOf(RULE) >= 0) {
    for (const [name, cut] of [['rows of a gated parent kept', 'function trioShut(s) { return gated(s); }'],
      ['rows of any child kept', 'function trioShut(s) { return false; }']]) {
      const res = await closedChecks(CORE.replace(RULE, cut), 'control');
      check('control: ' + name + ' fails the closed checks', res.some((c) => !c[1]));
    }
  }
  console.log(fails ? fails + ' check(s) failed' : 'All checks passed.');
  process.exit(fails ? 1 : 0);
})().catch((e) => { console.log('FAIL ' + (e && e.stack || e)); process.exit(1); });
