// map.js: Figure 1 of a person's page, the genome map. It draws the map document the server builds with
// `og map build` (GET map.json beside the page: gzip-wrapped columnar JSON, rust/crates/og-map/README.md):
// the 24 chromosomes as rows with their bands, the masked spans hatched, each finding with a place as its
// level's glyph, and the runs of homozygosity; a click opens a chromosome in place, where the bands, the
// findings and the runs get lanes of their own under a ruler. Every mark comes from the document, never from
// results.json. The margin's findings count, chromosome sentences, runs sentence and finding phrases are the
// document's (text, phrase, empty, note); this file adds only short lines of its own: an empty lane or
// stretch, the runs in view, how to move, and why a layer is not on the map.
//
// The veil. The document is the same for every reader: each row a closed result could reveal carries the id
// of that result (veil). Whether it is closed is the page's question (person.js veil()): a closed row is not
// drawn, has no hit target, no tooltip, no line in the margin and never reaches the address. The closed
// results are counted in a dashed chip that names no gene and no place, with an Open button that asks once,
// then opens the result through the page's own function, so the map and the card open together. A reload
// closes them again: the open state lives in the page's memory only.
//
// No library and no build step: one canvas, the camera {ci, c, ls, k} (open chromosome, centre bp, natural
// log of the span, 0 genome .. 1 open) eased toward its target each frame, one level() by bp per CSS pixel.
// Text only through textContent: no innerHTML, no style attribute. Everything hangs on window.PGPMap.
(function () {
  'use strict';
  var U = window.PGPUI;
  if (!U) { return; }
  var el = U.el, add = U.add;
  var LV = { act: '#d03b3b', medicine: '#ec835a', review: '#fab219', carrier: '#9085e9', trust: '#c3c2b7', info: '#c3c2b7' };
  var PRI = { act: 0, medicine: 1, review: 2, carrier: 3, trust: 4, info: 5 };
  var STAIN = { gneg: '#26272c', gpos25: '#323337', gpos50: '#3e3f42', gpos75: '#4c4c4d', gpos100: '#5a5955',
    gvar: '#25262a', stalk: '#25262a', acen: '#2e2f34' };
  var BRASS = '#c8ab78', INK = '#f0efec', INK2 = '#c3c2b7', INK3 = '#8f8d86', PAGE = '#0e0f10';
  var SANS = 'system-ui, -apple-system, "Segoe UI", sans-serif';
  var MONO = 'ui-monospace, "SF Mono", Menlo, monospace';
  // The levels, by bp per CSS pixel of the plot: chromosome above 16,667, region down to 833, then gene,
  // and reads below a 50 kb span (the read viewer's window is in bp, not pixels). The view never gets
  // narrower than 100 bp.
  var CHROM_BPP = 16667, REGION_BPP = 833, READS_SPAN = 50000, MIN_SPAN = 100;
  // The gauge's Region and Gene stops, in bp per CSS pixel, so each lands inside its level at any width.
  var REGION_STOP = 8000, GENE_STOP = 300;
  // A genome row's least height: 22 CSS px, 10 on a phone (the page there is long enough already).
  var ROW_MIN = 22, ROW_PHONE = 10;
  var reduced = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  function lerp(a, b, t) { return a + (b - a) * t; }
  function clamp(v, a, b) { return Math.max(a, Math.min(b, v)); }
  function mb(bp, d) { return (bp / 1e6).toFixed(d === undefined ? 1 : d); }
  function comma(n) { return String(Math.round(n)).replace(/\B(?=(\d{3})+(?!\d))/g, ','); }
  // A sentence of the document: pieces of text, the bold ones in <b>.
  function runs(list) {
    var p = el('p', {});
    (list || []).forEach(function (x) { add(p, x.b ? el('b', {}, x.t) : document.createTextNode(x.t)); });
    return p;
  }
  function short(n) { return String(n).replace(/^chr/, ''); }
  // A layer's states that put nothing on the map (the document's layer entry); any other state has landed.
  var NOT_ON_MAP = { missing: 1, unreadable: 1, not_run: 1, not_in_this_run: 1, running: 1, waiting: 1, failed: 1, skipped: 1 };
  function landed(lay) { return !!lay && !NOT_ON_MAP[lay.state]; }

  // --- what is drawn one by one: the one rule ---------------------------------------------------------
  // drawsAt(layer, row, level): whether the map draws this row of this layer as a mark of its own at this
  // level. The canvas, the whole-genome strip and the In view table all ask it, and nothing else in this file
  // decides what to draw. og_map::predicate (rust/crates/og-map/src/predicate.rs) holds the same rule, and the
  // browser check compares the two row by row (og map build --predicate-table). A finding draws at every level;
  // a run of homozygosity of 5 Mb or more at every level, a shorter one once its chromosome is open (at the
  // genome level it is part of the 1 Mb coverage bins); a gene from the region level down. The veil is apart:
  // a closed row is dropped before this rule is asked.
  var LEVELS = ['genome', 'chromosome', 'region', 'gene', 'reads'];
  function drawsAt(layer, row, level) {
    var i = LEVELS.indexOf(level);
    if (i < 0) { return false; }
    if (layer === 'findings') { return true; }
    if (layer === 'roh') { return i > 0 || !!(row && row.long); }
    if (layer === 'genes') { return i >= 2; }
    return false;
  }
  // A lane with more marks than this in view draws its density instead, and the table gives it one line.
  var DENSE = 4000;
  // A mark's key: the same for the canvas, the strip and the table, so a check can hold them to one set.
  function mark(layer, chrom, s, e, id) { return layer + '\t' + chrom + '\t' + s + '\t' + e + '\t' + (id || ''); }

  // --- the document: gzip -> JSON -> typed columns ------------------------------------------------------
  // Starts are gaps from the previous row (ds) and ends are lengths (l): s, e and the running maximum of e
  // are rebuilt here in one pass, so a view [a, b) is a binary search on that maximum.
  function columns(o) {
    var n = (o.ds || []).length, s = new Float64Array(n), e = new Float64Array(n), me = new Float64Array(n);
    var at = 0, m = 0;
    for (var i = 0; i < n; i++) {
      at += o.ds[i];
      s[i] = at;
      e[i] = at + o.l[i];
      m = Math.max(m, e[i]);
      me[i] = m;
    }
    return { n: n, s: s, e: e, me: me, o: o };
  }
  // The first row whose running maximum end passes a: every row from there with s < b meets [a, b).
  function first(c, a) {
    var lo = 0, hi = c.n;
    while (lo < hi) { var mid = (lo + hi) >> 1; if (c.me[mid] > a) { hi = mid; } else { lo = mid + 1; } }
    return lo;
  }
  function prepare(d, who) {
    var person = (d.people || [])[0] || {};
    var id = who || person.id;
    var per = (d.per || {})[id] || {};
    var chroms = ((d.ref && d.ref.chroms) || []).filter(function (c) { return !(c.n === 'chrY' && person.has_y === false); });
    var by = {};
    chroms.forEach(function (c, i) { by[c.n] = i; });
    var bands = {}, roh = {}, bins = {};
    chroms.forEach(function (c) {
      var b = d.ref.bands && d.ref.bands[c.n];
      if (b) { bands[c.n] = columns(b); }
      if (per.roh && per.roh[c.n]) { roh[c.n] = columns(per.roh[c.n]); }
      if (per.roh_bins && per.roh_bins[c.n]) { bins[c.n] = new Float64Array(per.roh_bins[c.n]); }
    });
    var f = per.findings || {}, veilIds = (d.dict && d.dict.veil) || [];
    var rows = [];
    (f.id || []).forEach(function (fid, i) {
      if (!(f.chrom[i] in by)) { return; }
      rows.push({ id: fid, level: f.level[i], gate: f.gate[i], chrom: f.chrom[i], pos: f.pos[i], s: f.s[i], e: f.e[i],
        gene: f.gene[i], title: f.title[i], href: f.href[i], phrase: f.phrase[i], veil: f.veil[i] >= 0 ? veilIds[f.veil[i]] : null,
        hgvs: (f.hgvs || [])[i] || '', zyg: (f.zyg || [])[i] || '', steps: (f.steps || [])[i] || [] });
    });
    var gated = {};
    rows.forEach(function (r) { if (r.gate) { gated[r.id] = r; } });
    // What can be closed: each veil tag, and each result to act on that carries another result's tag (a row
    // over two of them carries only the first one's tag, and one over APOE carries APOE's): its own id too.
    var closable = veilIds.slice();
    rows.forEach(function (r) { if (r.veil && r.level === 'act' && !r.gate && closable.indexOf(r.id) < 0) { closable.push(r.id); } });
    return { doc: d, id: id, label: person.label || id, chroms: chroms, by: by, bands: bands, roh: roh, bins: bins,
      binSize: (per.roh_bins && per.roh_bins.size) || 1e6, findings: rows, gated: gated, veilIds: veilIds, closable: closable,
      layers: d.layers || {}, text: d.text || {}, art: (d.dict && d.dict.art) || [], gencode: (d.ref && d.ref.gencode) || null };
  }

  // --- the figure ---------------------------------------------------------------------------------------
  // create({url, veil(id, gated) -> 'open'|'closed'|'opened', now(id) -> shown at this moment, open(id),
  // close(id), known(id) -> a finding of the page's results, closed() -> the ids of the results to act on the
  // page keeps closed at this moment, stepDone(step) -> ended at this moment, fallback() -> rows}): one figure
  // per page. The page calls refresh(key) on each of its redraws; a new key
  // fetches the document again and keeps the camera.
  function create(opts) {
    var M = null, key = null, failure = null, t0 = 0;
    var fig = el('figure', { 'class': 'row fig1 gmap', id: 'genome-map', 'aria-labelledby': 'gm-title' });
    var main = el('div', { 'class': 'gm-main' });
    var head = el('div', { 'class': 'gm-head' });
    var kicker = el('div', { 'class': 'fig-label' }, 'Figure 1 · Genome map');
    var title = el('h2', { 'class': 'gm-title', id: 'gm-title' }, 'The genome map');
    var gauge = el('div', { 'class': 'gm-gauge', role: 'group', 'aria-label': 'Zoom level' });
    var stopG = el('button', { type: 'button', 'data-stop': 'genome' }, 'Genome');
    var stopC = el('button', { type: 'button', 'data-stop': 'chromosome' }, 'Chromosome');
    var stopR = el('button', { type: 'button', 'data-stop': 'region' }, 'Region');
    var stopN = el('button', { type: 'button', 'data-stop': 'gene' }, 'Gene');
    add(gauge, stopG, stopC, stopR, stopN);
    // Zoom buttons for an open chromosome: the touch path to what the wheel and the + and - keys do.
    var zoom = el('div', { 'class': 'gm-zoom', role: 'group', 'aria-label': 'Zoom', hidden: 'hidden' });
    var zoomOut = el('button', { type: 'button', 'aria-label': 'Zoom out', 'data-zoom': 'out' }, '−');
    var zoomIn = el('button', { type: 'button', 'aria-label': 'Zoom in', 'data-zoom': 'in' }, '+');
    add(zoom, zoomOut, zoomIn);
    var crumb = el('nav', { 'class': 'gm-crumb mono', 'aria-label': 'Where the map is' });
    var chips = el('div', { 'class': 'gm-chips' });
    // Search (the / key): gene names, band names, the findings drawn and a place typed as the address takes it.
    var finder = el('div', { 'class': 'gm-find', role: 'search' });
    var findBox = el('input', { type: 'search', 'class': 'gm-find-box', id: 'gm-find', autocomplete: 'off', spellcheck: 'false',
      role: 'combobox', 'aria-expanded': 'false', 'aria-controls': 'gm-find-list', 'aria-autocomplete': 'list',
      'aria-label': 'Find a gene, a band, a finding or a place on the map', placeholder: 'Find a gene, band or finding  /' });
    var findList = el('ul', { 'class': 'gm-find-list', id: 'gm-find-list', role: 'listbox', 'aria-label': 'Places found', hidden: 'hidden' });
    add(finder, findBox, findList);
    add(head, kicker, add(el('div', { 'class': 'gm-titlerow' }), title, gauge, zoom), crumb, finder, chips);
    var stage = el('div', { 'class': 'gm-stage' });
    var canvas = el('canvas', { 'class': 'gm-canvas', tabindex: '0', role: 'img', 'aria-describedby': 'gm-narration' });
    var tip = el('div', { 'class': 'gm-tip', role: 'tooltip', hidden: 'hidden' });
    var toast = el('div', { 'class': 'gm-toast', role: 'status', hidden: 'hidden' });
    add(stage, canvas, tip, toast);
    // The marks drawn, as links: the keyboard and screen-reader path to each finding (Tab steps through them).
    var marks = el('ul', { 'class': 'gm-marks', 'aria-label': 'Findings on the map' });
    var fallback = el('div', { 'class': 'gm-fallback', hidden: 'hidden' });
    // The whole-genome strip over an open chromosome, and the In view table under the map (strip(), table()).
    var strip = el('canvas', { 'class': 'gm-strip', role: 'img', hidden: 'hidden' });
    var inview = el('details', { 'class': 'gm-inview' });
    add(main, head, strip, stage, inview, marks, fallback);
    var margin = el('aside', { 'class': 'gm-margin' });
    var mkick = el('div', { 'class': 'fig-label' }, 'The whole genome');
    var mline = el('div', { 'class': 'gm-mline mono' });
    var narr = el('div', { 'class': 'gm-narration', id: 'gm-narration', 'aria-live': 'polite' });
    // The card of the finding picked on the map, and the way to the nearest finding when none is in view.
    var card = el('section', { 'class': 'gm-card', 'aria-label': 'The finding picked on the map', hidden: 'hidden' });
    var nearBtn = el('button', { type: 'button', 'class': 'gm-nearest', hidden: 'hidden' });
    var paintBox = el('div', { 'class': 'gm-paint', role: 'group', 'aria-label': 'Paint the copies' });
    var layersBox = el('dl', { 'class': 'gm-layers' });
    var keyBox = el('div', { 'class': 'gm-key', 'aria-label': 'Key' });
    // The layers and the key: open beside the map, folded under it on a phone.
    var more = el('details', { 'class': 'gm-more' });
    more.open = !(window.matchMedia && window.matchMedia('(max-width: 700px)').matches);
    add(more, el('summary', {}, 'Layers and key'), paintBox, layersBox, keyBox);
    add(margin, mkick, mline, card, narr, nearBtn, more);
    add(fig, main, margin);

    var ctx = canvas.getContext('2d');
    var CW = 0, CH = 0, DPR = 1, phone = false;
    var cam = { ci: 0, c: 0, ls: 0, k: 0 }, tgt = { ci: 0, c: 0, ls: 0, k: 0 };
    var hover = { row: -1, hit: null }, hits = [], raf = 0, settleT = 0, hashT = 0, confirming = null, touched = false;
    // The results with no place on the map that its chip opened: each gets an "Opened in its section" chip while open.
    var inSection = [];
    var cache = null, cacheKey = '';
    // What the last frame drew one by one (mark() keys), the genome rows' share kept with their cache: the
    // browser check holds the In view table to it (canvas.pgpDrawn).
    var painted = [], cachePainted = [];
    var hatch = null;
    // The finding whose card the margin shows (an id of drawn() only), and the mark the keyboard is on.
    var sel = null, focusId = null;
    // The gene files, fetched once each: per chromosome when a view first reaches the region level, the
    // names when the search box gets its first key. Reference data, the same for every reader.
    var genes = {}, names = null;

    // --- the veil: what this page keeps closed at this moment ---
    function stateOf(id) { return opts.veil(id, !!(M && M.gated[id])); }
    // A row shows when its own result is not closed and the page shows its finding at this moment (a replay
    // before the finding's step ended does not).
    // A row that is a finding of the page (opts.known) follows the page's own veil only: the page already
    // hides every finding it keeps closed (opts.now), and its card shows the rest. Hiding such a row for
    // another result's tag would drop a finding whose card the page shows, and so place that closed result
    // inside it. A row the page has no finding for (the map's APOE row) follows its tag too.
    function shownRow(r) {
      var own = M.closable.indexOf(r.id) < 0 || stateOf(r.id) !== 'closed';
      var tag = (opts.known && opts.known(r.id)) || !r.veil || stateOf(r.veil) !== 'closed';
      return own && tag && (!opts.now || opts.now(r.id));
    }
    // The runs of homozygosity, as this page shows them now: a replay before their step ended has none yet.
    function rohLayer() {
      var lay = M.layers.roh;
      if (lay && landed(lay) && lay.step && opts.stepDone && !opts.stepDone(lay.step)) {
        return { state: 'waiting', step: lay.step, step_name: lay.step_name };
      }
      return lay;
    }
    function drawn() { return M ? M.findings.filter(shownRow) : []; }

    // --- geometry ---
    function chrom() { return M.chroms[cam.ci]; }
    // The genome rows keep their height while a phone's canvas grows for an open chromosome (size()).
    function rowH() { var lo = phone ? ROW_PHONE : ROW_MIN; return M ? clamp((genomeH() - 18) / M.chroms.length, lo, 44) : lo; }
    function genomeH() { return 18 + (M ? M.chroms.length : 24) * (phone ? ROW_PHONE : ROW_MIN); }
    function gutG() { return phone ? 30 : 40; }
    function pad() { return phone ? 10 : 22; }
    function gut() { return lerp(gutG(), phone ? 64 : 128, cam.k); }
    function plotW() { return CW - gut() - pad(); }
    function minSpan(len) { return Math.min(len, MIN_SPAN); }
    function openW() { return Math.max(200, CW - (phone ? 64 : 128) - pad()); }
    function view() {
      var span = Math.exp(cam.ls);
      return [cam.c - span / 2, cam.c + span / 2, span];
    }
    function gx(bp) { return gutG() + bp * (CW - gutG() - pad()) / M.chroms[0].len; }
    function ox(bp) { var v = view(); return gut() + (bp - v[0]) / v[2] * plotW(); }
    function xAt(bp, open) { return open ? lerp(gx(bp), ox(bp), cam.k) : gx(bp); }
    function rowY(i) { return 9 + i * rowH() + rowH() / 2; }
    // One rule for every lane, the gauge, the title and the address: the level of the camera shown now.
    function level() { return levelOf(cam); }
    function levelOf(c) {
      if (c.k < 0.5) { return 'genome'; }
      var span = Math.exp(c.ls), bpp = span / openW();
      if (bpp > CHROM_BPP) { return 'chromosome'; }
      if (bpp > REGION_BPP) { return 'region'; }
      return span < READS_SPAN ? 'reads' : 'gene';
    }
    function close_(lv) { return lv === 'gene' || lv === 'reads'; }
    // The lanes of an open chromosome, top to bottom: y of each, stretched to the canvas.
    function lanes() {
      var free = Math.max(0, CH - 230), st = 1 + Math.min(0.75, free / 400);
      var fy = 118, fh = 50 * Math.min(1.2, st), ry = fy + fh + 14, rh = 22 * st;
      var gy = ry + rh + 12, gh = clamp(CH - gy - 6, 60, 96);
      return { over: 10, ruler: 46, bands: 72, bh: 38, fy: fy, fh: fh, ry: ry, rh: rh, gy: gy, gh: gh };
    }

    // --- the canvas size, the hatch pattern ---
    // False while the stage has no width (not in the page yet, or hidden): sized then, a phone would get the
    // desktop's rows. The ResizeObserver draws once it has one.
    function size() {
      if (!stage.clientWidth) { return false; }
      var w = Math.max(280, Math.floor(stage.clientWidth));
      phone = w < 520;
      var n = M ? M.chroms.length : 24;
      // A phone's genome rows are short: an open chromosome makes the canvas taller, for its genes lane.
      var h = Math.max(18 + n * (phone ? ROW_PHONE : ROW_MIN), phone && tgt.k > 0 ? 330 : 0);
      DPR = Math.min(2, window.devicePixelRatio || 1);
      if (w !== CW || h !== CH || canvas.width !== Math.round(w * DPR)) {
        CW = w; CH = h;
        canvas.width = Math.round(w * DPR);
        canvas.height = Math.round(h * DPR);
        canvas.style.setProperty('width', w + 'px');
        canvas.style.setProperty('height', h + 'px');
        cacheKey = '';
      }
      return true;
    }
    function pattern(alpha) {
      var p = document.createElement('canvas');
      p.width = p.height = 6;
      var c = p.getContext('2d');
      c.strokeStyle = 'rgba(255,255,255,' + alpha + ')';
      c.lineWidth = 1;
      c.beginPath(); c.moveTo(-1, 7); c.lineTo(7, -1); c.moveTo(5, 7); c.lineTo(7, 5); c.moveTo(-1, 1); c.lineTo(1, -1); c.stroke();
      return ctx.createPattern(p, 'repeat');
    }

    // --- drawing ---
    function rounded(c, x, y, w, h, r) {
      r = Math.max(0, Math.min(r, w / 2, h / 2));
      c.moveTo(x + r, y); c.lineTo(x + w - r, y); c.arcTo(x + w, y, x + w, y + r, r); c.lineTo(x + w, y + h - r);
      c.arcTo(x + w, y + h, x + w - r, y + h, r); c.lineTo(x + r, y + h); c.arcTo(x, y + h, x, y + h - r, r);
      c.lineTo(x, y + r); c.arcTo(x, y, x + r, y, r);
    }
    // One chromosome's body: two rounded arms pinched at the centromere, filled with its bands, the masked
    // spans hatched. X maps a bp to a canvas x.
    function body(c, ch, X, yc, h, alpha) {
      var cen = ch.cen ? (ch.cen[0] + ch.cen[1]) / 2 : null;
      c.save();
      c.globalAlpha = alpha;
      c.beginPath();
      var x0 = X(0), x1 = X(ch.len), xm = cen === null ? null : X(cen);
      if (xm === null) { rounded(c, x0, yc - h / 2, x1 - x0, h, h / 2); } else {
        rounded(c, x0, yc - h / 2, Math.max(1, xm - x0 - 0.6), h, h / 2);
        rounded(c, xm + 0.6, yc - h / 2, Math.max(1, x1 - xm - 0.6), h, h / 2);
      }
      c.clip();
      c.fillStyle = STAIN.gneg;
      c.fillRect(x0, yc - h / 2, x1 - x0, h);
      var b = M.bands[ch.n];
      if (b) {
        for (var i = 0; i < b.n; i++) {
          var st = b.o.stain[i];
          if (st === 'gneg') { continue; }
          var a = X(b.s[i]), z = X(b.e[i]);
          if (z < 0 || a > CW) { continue; }
          c.fillStyle = STAIN[st] || STAIN.gneg;
          c.fillRect(a, yc - h / 2, Math.max(0.5, z - a), h);
        }
      }
      c.fillStyle = hatch;
      (ch.masked || []).forEach(function (m) { var a = X(m[0]), z = X(m[1]); c.fillRect(a, yc - h / 2, Math.max(1, z - a), h); });
      c.restore();
    }
    function glyph(c, lv, x, y, r) {
      c.save();
      c.fillStyle = LV[lv] || LV.info;
      c.strokeStyle = PAGE;
      c.lineWidth = 3;
      c.beginPath();
      if (lv === 'act') { c.moveTo(x, y - r); c.lineTo(x + r, y); c.lineTo(x, y + r); c.lineTo(x - r, y); c.closePath(); c.stroke(); c.fill(); }
      else if (lv === 'medicine') { rounded(c, x - r * 1.3, y - r * 0.7, r * 2.6, r * 1.4, r * 0.7); c.stroke(); c.fill(); }
      else if (lv === 'review') { c.moveTo(x, y - r); c.lineTo(x + r, y + r * 0.8); c.lineTo(x - r, y + r * 0.8); c.closePath(); c.stroke(); c.fill(); }
      else {
        var two = lv === 'carrier', rr = two ? r * 0.75 : r * 0.9;
        var at = two ? [x - rr * 0.75, x + rr * 0.75] : [x];
        at.forEach(function (cx) { c.moveTo(cx + rr, y); c.arc(cx, y, rr, 0, Math.PI * 2); });
        c.stroke();
        c.strokeStyle = LV[lv] || LV.info;
        c.lineWidth = 1.5;
        c.beginPath();
        at.forEach(function (cx) { c.moveTo(cx + rr, y); c.arc(cx, y, rr, 0, Math.PI * 2); });
        c.stroke();
      }
      c.restore();
    }
    function label(c, text, x, y, font, colour) {
      c.font = font;
      c.lineJoin = 'round';
      c.strokeStyle = PAGE;
      c.lineWidth = 4;
      c.strokeText(text, x, y);
      c.fillStyle = colour;
      c.fillText(text, x, y);
    }
    var widths = {};
    function measure(font, text) {
      var k = font + '\u0000' + text;
      if (!(k in widths)) { ctx.font = font; widths[k] = ctx.measureText(text).width; }
      return widths[k];
    }

    // The genome rows, drawn once into an offscreen canvas until the data, the open results or the size change.
    function genomeRows() {
      var k = [CW, CH, DPR, drawn().map(function (r) { return r.id; }).join(',')].join('|');
      if (cache && k === cacheKey) { return cache; }
      cache = cache || document.createElement('canvas');
      cache.width = canvas.width;
      cache.height = canvas.height;
      var c = cache.getContext('2d');
      c.setTransform(DPR, 0, 0, DPR, 0, 0);
      c.clearRect(0, 0, CW, CH);
      var fs = drawn().filter(function (r) { return drawsAt('findings', r, 'genome'); });
      var occupied = {};
      cachePainted = fs.map(function (r) { return mark('findings', r.chrom, r.s, r.e, r.id); });
      M.chroms.forEach(function (ch, i) {
        var yc = rowY(i);
        c.fillStyle = INK3;
        c.font = (phone ? '9px ' : '11px ') + SANS;
        c.textAlign = 'right';
        c.textBaseline = 'middle';
        c.fillText(short(ch.n), gutG() - 10, yc);
        c.textAlign = 'left';
        body(c, ch, gx, yc, 8, 1);
        // Runs of homozygosity under the body: the long ones solid, the rest as 1 Mb coverage bins.
        var rohOn = landed(rohLayer()), bins = rohOn && M.bins[ch.n], ry = yc + 6.5;
        if (bins) {
          c.fillStyle = INK2;
          for (var b = 0; b < bins.length; b++) {
            if (!bins[b]) { continue; }
            c.globalAlpha = 0.22 * Math.min(1, bins[b] / M.binSize);
            c.fillRect(gx(b * M.binSize), ry - 1, Math.max(1, gx((b + 1) * M.binSize) - gx(b * M.binSize)), 2);
          }
          c.globalAlpha = 1;
        }
        var rr = rohOn && M.roh[ch.n];
        if (rr) {
          c.strokeStyle = INK;
          c.lineWidth = 2.5;
          for (var j = 0; j < rr.n; j++) {
            if (!drawsAt('roh', { long: rr.o.long[j] }, 'genome')) { continue; }
            c.beginPath(); c.moveTo(gx(rr.s[j]), ry); c.lineTo(gx(rr.e[j]), ry); c.stroke();
            cachePainted.push(mark('roh', ch.n, rr.s[j], rr.e[j], ''));
          }
        }
        occupied[i] = [];
      });
      // Findings: glyphs pinned to their row, names to the right by level (act first), never over another
      // glyph or name of the row: a name that does not fit is left to the margin and the list.
      var byLevel = fs.slice().sort(function (a, b) { return PRI[a.level] - PRI[b.level] || a.s - b.s; });
      byLevel.forEach(function (r) {
        var i = M.by[r.chrom], yc = rowY(i), x = gx((r.s + r.e) / 2);
        if (gx(r.e) - gx(r.s) > 3) {
          c.fillStyle = LV[r.level] || LV.info;
          c.fillRect(gx(r.s), yc + 4.5, gx(r.e) - gx(r.s), 2);
        }
        glyph(c, r.level, x, yc - 1, phone ? 3.5 : 4.5);
        occupied[i].push([x - 6, x + 6]);
      });
      byLevel.forEach(function (r) {
        var i = M.by[r.chrom], yc = rowY(i), x = gx((r.s + r.e) / 2);
        var name = r.gene || '';
        if (!name) { return; }
        var font = '500 10.5px ' + SANS, w = measure(font, name), lx = x + 8;
        var clash = occupied[i].some(function (o) { return !(lx + w < o[0] || lx > o[1]) && !(o[0] === x - 6 && o[1] === x + 6); });
        if (lx + w > CW - 2 || clash) { return; }
        occupied[i].push([lx - 2, lx + w + 4]);
        c.textBaseline = 'middle';
        label(c, name, lx, yc - 1, font, INK2);
      });
      cacheKey = k;
      return cache;
    }

    function draw() {
      raf = 0;
      if (!M || failure) { return; }  // a failed reload leaves the old document undrawn
      var t = performance.now();
      if (!size()) { return; }  // the page may have changed width since the last frame (a phone turned, a fold opened)
      step();
      ctx.setTransform(1, 0, 0, 1, 0, 0);
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
      hits = [];
      painted = [];
      if (cam.k < 0.999) { drawGenome(); }
      if (cam.k > 0.001) { drawOpen(); }
      canvas.pgpDrawn = painted;
      drawStrip();
      if (!fig.hasAttribute('data-pgp-map-paint-ms') && t0) {
        fig.setAttribute('data-pgp-map-paint-ms', String(Math.round(performance.now() - t0)));
      }
      if (moving()) { schedule(); } else { settle(); }
      return t;
    }
    function drawGenome() {
      var k = cam.k, n = M.chroms.length;
      if (k < 0.001) {
        if (hover.row >= 0 && level() === 'genome') {
          ctx.fillStyle = 'rgba(255,255,255,.035)';
          ctx.fillRect(0, rowY(hover.row) - rowH() / 2, CW, rowH());
          var ch = M.chroms[hover.row], words = 'Open chromosome ' + short(ch.n) + ' ›';
          var font = '500 11px ' + SANS, w = measure(font, words), x = gx(ch.len) + 12;
          if (x + w < CW - 4) { ctx.textBaseline = 'middle'; label(ctx, words, x, rowY(hover.row), font, BRASS); }
        }
        ctx.drawImage(genomeRows(), 0, 0, CW, CH);
        painted = painted.concat(cachePainted);
        M.chroms.forEach(function (c, i) {
          hits.push({ kind: 'row', i: i, x0: 0, x1: CW, y0: rowY(i) - rowH() / 2, y1: rowY(i) + rowH() / 2 });
        });
        drawn().forEach(function (r) {
          var x = gx((r.s + r.e) / 2), y = rowY(M.by[r.chrom]) - 1;
          if (r.id === sel || r.id === focusId) {
            ctx.strokeStyle = BRASS;
            ctx.lineWidth = 2;
            ctx.beginPath(); ctx.arc(x, y, 8, 0, Math.PI * 2); ctx.stroke();
          }
          hits.push({ kind: 'finding', r: r, x: x, y: y, rad: 7 });
        });
        return;
      }
      // The motion: the other rows slide away from the one opening and fade.
      var a = Math.pow(1 - k, 1.4);
      for (var i = 0; i < n; i++) {
        if (i === cam.ci) { continue; }
        var dy = k * CH * 0.9 * (0.6 + 0.04 * Math.abs(i - cam.ci)) * (i < cam.ci ? -1 : 1);
        ctx.save();
        ctx.globalAlpha = a;
        ctx.translate(0, dy);
        ctx.fillStyle = INK3;
        ctx.font = '11px ' + SANS;
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        ctx.fillText(short(M.chroms[i].n), gutG() - 10, rowY(i));
        ctx.textAlign = 'left';
        body(ctx, M.chroms[i], gx, rowY(i), 8, 1);
        ctx.restore();
      }
    }
    // A place on the ruler, in Mb with as many decimals as the step needs, in bp below a 1 kb step.
    function tick(at, stp) {
      if (stp < 1000) { return comma(at); }
      return mb(at, stp >= 1e6 ? 0 : Math.ceil(-Math.log10(stp / 1e6) - 1e-9)) + ' Mb';
    }
    // The reading line's place: closer views say it more precisely.
    function where(bp, span) {
      if (span < READS_SPAN) { return comma(bp); }
      return mb(bp, span >= 20e6 ? 1 : span >= 2e6 ? 2 : 3) + ' Mb';
    }
    function niceStep(span, px) {
      var raw = span / Math.max(1, px / 90), p = Math.pow(10, Math.floor(Math.log10(raw)));
      var m = raw / p;
      return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
    }
    function drawOpen() {
      var k = cam.k, ch = chrom(), L = lanes(), v = view(), g = gut(), W = plotW();
      var X = function (bp) { return xAt(bp, true); };
      var yb = lerp(rowY(cam.ci), L.bands + L.bh / 2, k), hb = lerp(8, L.bh, k);
      ctx.save();
      // The open row's label moves into the gutter.
      ctx.fillStyle = INK3;
      ctx.font = '11px ' + SANS;
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      ctx.fillText(k > 0.5 ? 'Bands' : short(ch.n), lerp(gutG() - 10, g - 12, k), yb);
      ctx.restore();
      ctx.save();
      ctx.beginPath();
      ctx.rect(lerp(0, g - 2, k), 0, CW, CH);
      ctx.clip();
      body(ctx, ch, X, yb, hb, 1);
      ctx.restore();
      ctx.save();
      ctx.globalAlpha = k;
      // Overview: the whole chromosome, the view as a brass window, a funnel to the ruler.
      var ow = W, oX = function (bp) { return g + bp / ch.len * ow; };
      body(ctx, ch, oX, L.over + 6, 8, 0.8);
      ctx.strokeStyle = BRASS;
      ctx.lineWidth = 1.5;
      var wa = oX(clamp(v[0], 0, ch.len)), wz = oX(clamp(v[1], 0, ch.len));
      ctx.strokeRect(wa, L.over, Math.max(3, wz - wa), 12);
      ctx.globalAlpha = k * 0.35;
      ctx.beginPath(); ctx.moveTo(wa, L.over + 12); ctx.lineTo(g, L.ruler); ctx.moveTo(wz, L.over + 12); ctx.lineTo(g + W, L.ruler); ctx.stroke();
      ctx.globalAlpha = k;
      // The ruler, in Mb.
      ctx.strokeStyle = 'rgba(255,255,255,.14)';
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(g, L.ruler); ctx.lineTo(g + W, L.ruler); ctx.stroke();
      var stp = niceStep(v[2], W);
      ctx.fillStyle = INK3;
      ctx.font = '10.5px ' + SANS;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      for (var at = Math.ceil(v[0] / stp) * stp; at <= v[1]; at += stp) {
        var x = ox(at);
        ctx.beginPath(); ctx.moveTo(x, L.ruler); ctx.lineTo(x, L.ruler + 5); ctx.stroke();
        ctx.fillText(tick(at, stp), x, L.ruler + 8);
      }
      // Band names where they fit: the arm's name at the chromosome level, with its chromosome closer in.
      var b = M.bands[ch.n], lvl = level();
      if (b) {
        ctx.font = '10px ' + SANS;
        ctx.textBaseline = 'middle';
        for (var i = first(b, v[0]); i < b.n && b.s[i] < v[1]; i++) {
          var a = ox(b.s[i]), z = ox(b.e[i]), nm = b.o.name[i] && lvl !== 'chromosome' ? short(ch.n) + b.o.name[i] : b.o.name[i];
          if (nm && z - a > measure('10px ' + SANS, nm) + 6) {
            ctx.fillStyle = ['gpos75', 'gpos100'].indexOf(b.o.stain[i]) >= 0 ? INK : INK2;
            ctx.fillText(nm, (Math.max(a, g) + Math.min(z, g + W)) / 2, L.bands + L.bh / 2);
          }
          hits.push({ kind: 'band', x0: a, x1: z, y0: L.bands, y1: L.bands + L.bh, name: nm, s: b.s[i], e: b.e[i] });
        }
      }
      (ch.masked || []).forEach(function (m) {
        hits.push({ kind: 'masked', x0: ox(m[0]), x1: ox(m[1]), y0: L.bands, y1: L.bands + L.bh, s: m[0], e: m[1] });
      });
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      ctx.fillStyle = INK3;
      ctx.font = '11px ' + SANS;
      ctx.fillText('Findings', g - 12, L.fy + 16);
      ctx.fillText(phone ? 'Runs' : 'Homozygous runs', g - 12, L.ry + L.rh / 2);
      ctx.fillText('Genes', g - 12, L.gy + 12);
      ctx.textAlign = 'left';
      ctx.save();
      ctx.beginPath(); ctx.rect(g - 6, 0, W + 12, CH); ctx.clip();
      openGenes(L, v, ch);
      openFindings(L, v, ch);
      openRoh(L, v, ch);
      // The reading line: a dashed brass line through every lane at the centre, its place in a pill.
      var cx = g + W / 2;
      ctx.setLineDash([3, 3]);
      ctx.strokeStyle = BRASS;
      ctx.globalAlpha = k * 0.8;
      ctx.beginPath(); ctx.moveTo(cx, L.ruler + 22); ctx.lineTo(cx, L.gy + L.gh); ctx.stroke();
      ctx.setLineDash([]);
      ctx.restore();
      var pill = where(cam.c, v[2]), pw = measure('600 10.5px ' + SANS, pill) + 12;
      ctx.fillStyle = BRASS;
      rounded(ctx, cx - pw / 2, L.ruler - 20, pw, 15, 7.5);
      ctx.beginPath(); rounded(ctx, cx - pw / 2, L.ruler - 20, pw, 15, 7.5); ctx.fill();
      ctx.fillStyle = PAGE;
      ctx.font = '600 10.5px ' + SANS;
      ctx.textAlign = 'center';
      ctx.fillText(pill, cx, L.ruler - 12.5);
      ctx.restore();
    }
    function laneWords(c, text, x, y) {
      c.fillStyle = INK3;
      c.font = 'italic 11.5px ' + SANS;
      c.textAlign = 'left';
      c.textBaseline = 'middle';
      c.fillText(text, x, y);
    }
    function openFindings(L, v, ch) {
      var lay = M.layers.findings || {}, g = gut();
      if (!landed(lay)) { laneWords(ctx, layerWords(lay), g + 4, L.fy + 16); return; }
      var fs = drawn().filter(function (r) { return r.chrom === ch.n; });
      var lvl = level();
      var inView = fs.filter(function (r) { return r.e > v[0] && r.s < v[1] && drawsAt('findings', r, lvl); });
      if (!inView.length) {
        var nr = nearest(), words0 = fs.length ? 'No finding in this stretch.' : noFindingHere();
        laneWords(ctx, words0, g + 4, L.fy + 16);
        if (nr) {
          // The way to the nearest finding, on the side it lies: the same button as the margin's.
          var t = nr.words, font = '500 11.5px ' + SANS, tw = measure(font, t), ww = measure('italic 11.5px ' + SANS, words0);
          var tx = nr.side < 0 ? g + 4 : Math.max(g + ww + 18, g + plotW() - tw - 4), ty = L.fy + 34;
          ctx.textBaseline = 'middle';
          label(ctx, t, tx, ty, font, BRASS);
          hits.push({ kind: 'nearest', x0: tx - 4, x1: tx + tw + 4, y0: ty - 9, y1: ty + 9 });
        }
        return;
      }
      var used = [], close = close_(level());
      inView.sort(function (a, b) { return PRI[a.level] - PRI[b.level] || a.s - b.s; }).forEach(function (r) {
        var x = ox((r.s + r.e) / 2), y = L.fy + 14;
        if (r.e - r.s > v[2] / 2000) {
          ctx.fillStyle = LV[r.level] || LV.info;
          ctx.fillRect(ox(r.s), y + 9, Math.max(2, ox(r.e) - ox(r.s)), 2);
        }
        ctx.strokeStyle = LV[r.level] || LV.info;
        ctx.globalAlpha = 0.5;
        ctx.lineWidth = 1.2;
        ctx.beginPath(); ctx.moveTo(x, y + 8); ctx.lineTo(x, L.bands + L.bh); ctx.stroke();
        if (close) {
          // Close in, the change is a dashed tick down through every lane, so its place in the gene shows.
          ctx.setLineDash([4, 3]);
          ctx.beginPath(); ctx.moveTo(x, y + 8); ctx.lineTo(x, L.gy + L.gh); ctx.stroke();
          ctx.setLineDash([]);
        }
        ctx.globalAlpha = cam.k;
        if (r.id === sel || r.id === focusId) {
          ctx.strokeStyle = BRASS;
          ctx.lineWidth = 2;
          ctx.beginPath(); ctx.arc(x, y, 10, 0, Math.PI * 2); ctx.stroke();
        }
        glyph(ctx, r.level, x, y, 6);
        hits.push({ kind: 'finding', r: r, x: x, y: y, rad: 9 });
        painted.push(mark('findings', r.chrom, r.s, r.e, r.id));
        // The level's words under the name; close in, the change itself (HGVS) in mono when it has one.
        var name = r.gene || r.title, mono = close && r.hgvs, words = mono ? r.hgvs : U.LEVEL_WORDS[r.level] || r.level;
        var wf = mono ? '10.5px ' + MONO : '10.5px ' + SANS;
        var w = Math.max(measure('500 12px ' + SANS, name), measure(wf, words));
        var lx = x + 11;
        if (lx + w > CW - 4 || used.some(function (u) { return !(lx + w < u[0] || lx > u[1]); })) { return; }
        used.push([lx - 4, lx + w + 6]);
        ctx.textBaseline = 'middle';
        label(ctx, name, lx, y - 1, '500 12px ' + SANS, INK);
        label(ctx, words, lx, y + 13, wf, INK3);
      });
    }
    function openRoh(L, v, ch) {
      var lay = rohLayer() || {}, g = gut(), y = L.ry + L.rh / 2;
      if (!landed(lay)) { laneWords(ctx, layerWords(lay), g + 4, y); return; }
      var r = M.roh[ch.n];
      if (!r || !r.n) { laneWords(ctx, 'No run of homozygosity on this chromosome.', g + 4, y); return; }
      var lvl = level();
      for (var i = first(r, v[0]); i < r.n && r.s[i] < v[1]; i++) {
        // A run that ends before the view (a short one after a longer one) is out of view: not drawn.
        if (r.e[i] <= v[0] || !drawsAt('roh', { long: r.o.long[i] }, lvl)) { continue; }
        painted.push(mark('roh', ch.n, r.s[i], r.e[i], ''));
        var a = ox(r.s[i]), z = ox(r.e[i]), art = M.art[r.o.art[i]] || '';
        if (r.o.long[i]) {
          ctx.fillStyle = INK2;
          ctx.fillRect(a, y - 3, Math.max(1, z - a), 6);
          var t = mb(r.e[i] - r.s[i]) + ' Mb run of homozygosity', font = '11px ' + SANS;
          if (measure(font, t) + 8 > z - a) { t = mb(r.e[i] - r.s[i]) + ' Mb'; }
          if (measure(font, t) + 8 < z - a) { ctx.textBaseline = 'middle'; label(ctx, t, a + 4, y - 11, font, INK2); }
        } else if (art) {
          ctx.fillStyle = hatch;
          ctx.fillRect(a, y - 3, Math.max(1, z - a), 6);
        } else {
          ctx.fillStyle = 'rgba(195,194,183,.35)';
          ctx.fillRect(a, y - 0.5, Math.max(1, z - a), 1);
        }
        hits.push({ kind: 'roh', x0: a, x1: Math.max(a + 3, z), y0: y - 6, y1: y + 6, s: r.s[i], e: r.e[i], art: art, long: r.o.long[i] });
      }
    }
    // --- the genes lane: from the region level down, the chromosome's gene file (og map genes) ---
    // Beside map.json on a person's page: /runs/<run>/map.json, /map/genes/<chrom>.json.
    function genesUrl(name) {
      if (opts.genes) { return opts.genes + name + '.json'; }
      return new URL('../../map/genes/' + encodeURIComponent(name) + '.json', new URL(opts.url, window.location.href)).href;
    }
    // A gzip-wrapped JSON file, or an error carrying the server's first line (og's reason on a 404).
    function fetchGz(url) {
      return fetch(url, { credentials: 'same-origin' }).then(function (r) {
        if (!r.ok) {
          return r.text().then(function (t) { throw new Error((t || '').trim().split('\n')[0] || ('the server answered ' + r.status)); });
        }
        return new Response(r.body.pipeThrough(new window.DecompressionStream('gzip'))).text();
      }).then(function (t) { return JSON.parse(t); });
    }
    function genesOf(ch) {
      if (!M.gencode) { return { state: 'none' }; }
      var gs = genes[ch.n];
      if (gs) { return gs; }
      gs = genes[ch.n] = { state: 'loading' };
      fetchGz(genesUrl(ch.n)).then(function (d) {
        gs.g = columns(d.g || {});
        gs.x = columns(d.x || {});
        gs.state = 'ready';
        schedule();
        settle();
      }, function (err) {
        gs.state = 'failed';
        gs.note = String((err && err.message) || err).replace(/\.$/, '');
        schedule();
        settle();
      });
      return gs;
    }
    // The gene under BP: the gene level's title.
    function geneAt(ch, bp) {
      var gs = genes[ch.n];
      if (!gs || gs.state !== 'ready') { return null; }
      var G = gs.g, best = -1;
      // A protein-coding gene with a name first, then any protein-coding one, then any named one.
      var rank = function (j) { return (G.o.type[j] === 'protein_coding' ? 0 : 2) + (/^ENSG\d/.test(G.o.name[j]) ? 1 : 0); };
      for (var i = first(G, bp); i < G.n && G.s[i] <= bp; i++) {
        if (G.e[i] <= bp) { continue; }
        if (best < 0 || rank(i) < rank(best) || rank(i) === rank(best) && G.e[i] - G.s[i] < G.e[best] - G.s[best]) { best = i; }
      }
      return best < 0 ? null : G.o.name[best];
    }
    var occ = [], genesDrawn = 0;
    // One row of the lane is taken from a to z (CSS px)? Rows are pixel maps, so placing a gene costs its width.
    function taken(row, a, z) {
      for (var x = Math.max(0, Math.floor(a)), end = Math.min(CW - 1, Math.ceil(z)); x <= end; x++) { if (row[x]) { return true; } }
      return false;
    }
    function take(row, a, z) { row.fill(1, Math.max(0, Math.floor(a)), Math.min(CW, Math.ceil(z) + 1)); }
    function openGenes(L, v, ch) {
      var g = gut(), lv = level(), y0 = L.gy, ty = y0 + 12;
      genesDrawn = 0;
      if (!drawsAt('genes', null, lv)) { laneWords(ctx, 'Zoom in to a region to see its genes.', g + 4, ty); return; }
      var gs = genesOf(ch);
      if (gs.state === 'none') { laneWords(ctx, 'Not on the map: gene names need step 08’s GENCODE file.', g + 4, ty); return; }
      if (gs.state === 'loading') { laneWords(ctx, 'Loading the genes…', g + 4, ty); return; }
      if (gs.state === 'failed') { laneWords(ctx, 'Not on the map: ' + gs.note + '.', g + 4, ty); return; }
      var G = gs.g, X = gs.x, list = [];
      for (var i = first(G, v[0]); i < G.n && G.s[i] < v[1]; i++) { if (G.e[i] > v[0]) { list.push(i); } }
      if (!list.length) { laneWords(ctx, 'No gene in this stretch.', g + 4, ty); return; }
      var W = plotW();
      if (list.length > DENSE) {
        // Too many to draw one by one: how many genes start in each 2 px, and the words to zoom in.
        var bins = new Uint16Array(Math.ceil(W / 2) + 1), top = 1;
        list.forEach(function (j) { var b = Math.floor((ox(Math.max(G.s[j], v[0])) - g) / 2); if (b >= 0 && b < bins.length) { top = Math.max(top, ++bins[b]); } });
        ctx.fillStyle = INK2;
        for (var b = 0; b < bins.length; b++) {
          if (!bins[b]) { continue; }
          ctx.globalAlpha = 0.25 + 0.6 * bins[b] / top;
          ctx.fillRect(g + b * 2, y0 + 26, 2, 22);
        }
        ctx.globalAlpha = 1;
        laneWords(ctx, 'Dense here: zoom in to see each one.', g + 4, ty);
        return;
      }
      var close = close_(lv), RH = 29, fg = {}, labels = 0;
      drawn().forEach(function (r) { if (r.chrom === ch.n && r.gene) { fg[r.gene] = 1; } });
      // Who gets a place first: the genes a finding names, then the protein-coding ones, then the rest,
      // each group left to right. Three rows, each a strip for the bars and a strip for the names under
      // them: a gene takes the first row whose bar strip is free, and its name goes under it when the name
      // strip there is free too. A gene with no free bar strip is drawn over the last row, without a name.
      list.sort(function (a, b) {
        var ra = fg[G.o.name[a]] ? 0 : G.o.type[a] === 'protein_coding' ? 1 : 2, rb = fg[G.o.name[b]] ? 0 : G.o.type[b] === 'protein_coding' ? 1 : 2;
        return ra - rb || G.s[a] - G.s[b];
      });
      for (var r0 = 0; r0 < 6; r0++) {
        if (!occ[r0] || occ[r0].length !== CW) { occ[r0] = new Uint8Array(CW); } else { occ[r0].fill(0); }
      }
      list.forEach(function (j) {
        var name = G.o.name[j], bold = !!fg[name], font = (bold ? '600 ' : '') + '11px ' + SANS;
        var xa = Math.max(g - 4, ox(G.s[j])), xz = Math.min(g + W + 4, Math.max(xa + 1, ox(G.e[j])));
        var w = measure(font, name), lx = clamp(xa, g, g + W - w), row = -1;
        for (var k = 0; k < 3 && row < 0; k++) { if (!taken(occ[k], xa - 1, xz + 1)) { row = k; } }
        if (row < 0) { row = 2; } else { take(occ[row], xa - 1, xz + 1); }
        var named = labels < 80 && !taken(occ[row + 3], lx - 4, lx + w + 6);
        if (named) { take(occ[row + 3], lx - 4, lx + w + 6); }
        genesDrawn++;
        painted.push(mark('genes', ch.n, G.s[j], G.e[j], name));
        var yc = y0 + 9 + row * RH;
        ctx.fillStyle = ctx.strokeStyle = bold ? INK : G.o.type[j] === 'protein_coding' ? INK2 : INK3;
        if (close) {
          // A gene model: the span as a 1 px line with its strand's chevrons, the exon blocks inside it.
          ctx.fillRect(xa, yc - 0.5, xz - xa, 1);
          var dir = G.o.strand[j] === '+' ? 1 : G.o.strand[j] === '-' ? -1 : 0;
          if (dir) {
            ctx.lineWidth = 1;
            ctx.beginPath();
            for (var cx = xa + 14; cx < xz - 6; cx += 28) { ctx.moveTo(cx - 2 * dir, yc - 3); ctx.lineTo(cx + 2 * dir, yc); ctx.lineTo(cx - 2 * dir, yc + 3); }
            ctx.stroke();
          }
          for (var x = first(X, G.s[j]); x < X.n && X.s[x] < G.e[j]; x++) {
            var ea = Math.max(X.s[x], G.s[j]), ez = Math.min(X.e[x], G.e[j]);
            if (ez <= ea || ez <= v[0] || ea >= v[1]) { continue; }
            var a = ox(ea);
            ctx.fillRect(a, yc - 5, Math.max(1, ox(ez) - a), 10);
          }
        } else {
          ctx.fillRect(xa, yc - 2.5, xz - xa, 5);
        }
        if (named) {
          labels++;
          ctx.textBaseline = 'middle';
          label(ctx, name, lx, yc + 13, font, bold ? INK : INK2);
        }
        hits.push({ kind: 'gene', x0: xa, x1: Math.max(xz, xa + 3), y0: yc - 7, y1: yc + 18, name: name, s: G.s[j], e: G.e[j],
          strand: G.o.strand[j], type: G.o.type[j] });
      });
    }
    function noFindingHere() {
      return closedAny() ? 'No open finding with a place on this chromosome.' : 'No finding with a place on this chromosome.';
    }
    function closedAny() { return closedActs().length > 0 || M.closable.some(function (id) { return stateOf(id) === 'closed'; }); }
    // The results to act on this page keeps closed, in the order Open opens them: the page's own order (the
    // order of its closed cards), with a place on the map or without. Its length is what the chip counts, so
    // neither the count nor the order tells whether a closed result has a place. A page that passes no list
    // closes only the map's own rows, in map order.
    function closedActs() {
      var page = opts.closed ? opts.closed() : null;
      if (page) { return page.filter(function (id) { return !M.gated[id]; }); }
      return M.closable.filter(function (id) { return !M.gated[id] && stateOf(id) === 'closed'; });
    }
    function layerWords(lay) {
      if (lay.state === 'missing') { return 'Not on the map: ' + (lay.note || ('step ' + lay.step + ' wrote no ' + (lay.source || 'file') + '.')); }
      if (lay.state === 'unreadable') { return 'Not on the map: ' + (lay.source || '') + ' could not be read' + (lay.note ? ' (' + lay.note + ')' : '') + '.'; }
      if (lay.state === 'not_run') { return 'Not on the map: step ' + lay.step + ' was not part of this run.'; }
      if (lay.state === 'not_in_this_run') { return 'Not on the map: not part of this run.'; }
      if (lay.state === 'running' || lay.state === 'waiting') { return 'Not on the map yet: step ' + lay.step + ' has not ended.'; }
      if (lay.state === 'failed') { return 'Not on the map: step ' + lay.step + ' failed.'; }
      if (lay.state === 'skipped') { return 'Not on the map: step ' + lay.step + ' was skipped.'; }
      return lay.note || '';
    }

    // --- the camera ---
    function moving() {
      return Math.abs(cam.k - tgt.k) > 0.002 || Math.abs(cam.ls - tgt.ls) > 0.002 || Math.abs(cam.c - tgt.c) > Math.exp(cam.ls) * 0.0005;
    }
    var last = 0;
    function step() {
      var now = performance.now(), dt = last ? Math.min(64, now - last) : 16;
      last = now;
      if (reduced || cam.ci !== tgt.ci && cam.k < 0.01) { cam.ci = tgt.ci; }
      if (reduced) { cam.c = tgt.c; cam.ls = tgt.ls; cam.k = tgt.k; return; }
      if (cam.ci !== tgt.ci) {
        // Another chromosome: close the open one first, then open the new one.
        cam.k += (0 - cam.k) * (1 - Math.exp(-dt / 60));
        if (cam.k < 0.01) { cam.k = 0; cam.ci = tgt.ci; cam.c = tgt.c; cam.ls = tgt.ls; }
        return;
      }
      cam.c += (tgt.c - cam.c) * (1 - Math.exp(-dt / 75));
      cam.ls += (tgt.ls - cam.ls) * (1 - Math.exp(-dt / 75));
      cam.k += (tgt.k - cam.k) * (1 - Math.exp(-dt / 95));
      if (!moving()) { cam.c = tgt.c; cam.ls = tgt.ls; cam.k = tgt.k; }
    }
    function schedule() { if (!raf && !failure) { raf = window.requestAnimationFrame(function () { last = last || performance.now(); draw(); }); } }
    function clampTgt() {
      var ch = M.chroms[tgt.ci], span = clamp(Math.exp(tgt.ls), minSpan(ch.len), ch.len);
      tgt.ls = Math.log(span);
      tgt.c = clamp(tgt.c, span / 2, ch.len - span / 2);
    }
    function openChrom(i, c, span, push) {
      if (!M || i < 0 || i >= M.chroms.length) { return; }
      var ch = M.chroms[i];
      tgt.ci = i;
      tgt.k = 1;
      tgt.c = c === undefined ? ch.len / 2 : c;
      tgt.ls = Math.log(span === undefined ? ch.len : span);
      clampTgt();
      if (cam.k < 0.01) { cam.ci = i; cam.c = tgt.c; cam.ls = tgt.ls; }
      changed(push);
    }
    function closeChrom(push) {
      tgt.k = 0;
      changed(push);
    }
    function changed(push) {
      touched = true;
      last = 0;
      schedule();
      writeHash(push);
    }

    // --- the address: #genome, #chr7, #chr7:117.4-137.4M; never a closed result's id, gene or place ---
    function locusWords() {
      if (tgt.k < 0.5) { return 'genome'; }
      var ch = M.chroms[tgt.ci], span = Math.exp(tgt.ls);
      if (span >= ch.len * 0.999) { return ch.n; }
      var a = Math.max(0, tgt.c - span / 2), z = Math.min(ch.len, tgt.c + span / 2);
      if (span > 2e6) { return ch.n + ':' + mb(a) + '-' + mb(z) + 'M'; }
      if (span > 6e4) { return ch.n + ':' + mb(a, 3) + '-' + mb(z, 3) + 'M'; }
      return ch.n + ':' + Math.round(a + 1) + '-' + Math.round(z);
    }
    function writeHash(push) {
      if (!M || !window.history || !window.history.replaceState) { return; }
      window.clearTimeout(hashT);
      var go = function () {
        var h = '#' + locusWords();
        if (window.location.hash === h) { return; }
        if (push) { window.history.pushState(null, '', h); } else { window.history.replaceState(null, '', h); }
      };
      if (push) { go(); } else { hashT = window.setTimeout(go, 250); }
    }
    // A place in the address: 'genome', 'chr7', 'chr7:117.4-117.7M', 'chr7:120k-140k', 'chr7:117480026-117668665'
    // (1-based, inclusive), or the aliases view=genome, view=chromosome&chr=7, chr=7. Returns null when the
    // address is not a map place, false when it is one that cannot be read.
    function parsePlace(h) {
      var m, kv = {};
      if (/^genome$/i.test(h)) { return { k: 0 }; }
      if (/^(view|chr)=/.test(h)) {
        h.split('&').forEach(function (p) { var x = p.split('='); kv[x[0]] = decodeURIComponent(x[1] || ''); });
        if (kv.view === 'genome') { return { k: 0 }; }
        if (!kv.chr) { return kv.view ? false : null; }
        h = 'chr' + String(kv.chr).replace(/^chr/i, '');
      }
      if (!/^chr/i.test(h)) { return null; }
      var name = 'chr' + h.slice(3).split(':')[0].toUpperCase().replace(/^0+/, '');
      if (!(name in M.by)) { return false; }
      var ch = M.chroms[M.by[name]], rest = h.indexOf(':') > 0 ? h.slice(h.indexOf(':') + 1) : '';
      if (!rest) { return { k: 1, ci: M.by[name] }; }
      if ((m = /^([0-9.]+)-([0-9.]+)([MmKk])$/.exec(rest))) {
        var u = /m/i.test(m[3]) ? 1e6 : 1e3, a = Number(m[1]) * u, z = Number(m[2]) * u;
        if (!(z > a)) { return false; }
        return { k: 1, ci: M.by[name], s: a, e: z };
      }
      if ((m = /^([0-9,]+)-([0-9,]+)$/.exec(rest))) {
        var s = Number(m[1].replace(/,/g, '')) - 1, e = Number(m[2].replace(/,/g, ''));
        if (!(e > s)) { return false; }
        return { k: 1, ci: M.by[name], s: s, e: e };
      }
      void ch;
      return false;
    }
    function readHash(scroll) {
      if (!M) { return; }
      var h;
      try { h = decodeURIComponent((window.location.hash || '').slice(1)); } catch (e) { return; }
      // Back to the address with no place (the page as first loaded): the whole genome again.
      if (!h) { if (!scroll && tgt.k) { tgt.k = 0; last = 0; schedule(); } return; }
      var p = parsePlace(h);
      if (p === null) { return; }
      if (p === false) { say('Could not read the place “' + h + '”'); return; }
      if (p.k === 0) { tgt.k = 0; } else {
        var ch = M.chroms[p.ci];
        openChrom(p.ci, p.s === undefined ? ch.len / 2 : (p.s + p.e) / 2, p.s === undefined ? ch.len : p.e - p.s, false);
        window.clearTimeout(hashT);
      }
      if (scroll) {
        cam.ci = tgt.ci; cam.k = tgt.k; cam.c = tgt.c; cam.ls = tgt.ls;
        fig.scrollIntoView({ block: 'start' });
      }
      schedule();
    }
    function say(words) {
      toast.textContent = words;
      toast.hidden = false;
      window.clearTimeout(say.t);
      say.t = window.setTimeout(function () { toast.hidden = true; }, 4000);
    }

    // --- the margin: settled, not per frame ---
    function settle() {
      window.clearTimeout(settleT);
      // A gene file that lands while the camera moves asks for a settle too: the last frame of the motion settles
      // again, so the margin and the table describe the view the camera comes to rest on, not one it passes.
      settleT = window.setTimeout(function () { if (!moving()) { describe(); } }, 140);
    }
    function describe() {
      if (!M || failure) { return; }
      var lv = level(), fs = drawn(), ch = chrom();
      fig.setAttribute('data-level', lv);
      fig.setAttribute('data-chrom', lv === 'genome' ? '' : ch.n);
      fig.setAttribute('data-view', locusWords());
      fig.setAttribute('data-genes', lv === 'genome' || lv === 'chromosome' ? '' : String(genesDrawn));
      stopG.setAttribute('aria-pressed', String(lv === 'genome'));
      stopC.setAttribute('aria-pressed', String(lv === 'chromosome'));
      stopR.setAttribute('aria-pressed', String(lv === 'region'));
      stopN.setAttribute('aria-pressed', String(close_(lv)));
      zoom.hidden = lv === 'genome';
      crumb.textContent = '';
      crumb.hidden = lv === 'genome';  // the gauge already says Genome
      var gl = el('button', { type: 'button', 'class': 'gm-crumb-g' }, 'Genome');
      gl.addEventListener('click', function () { closeChrom(true); });
      add(crumb, gl);
      var nY = M.chroms.length;
      if (lv === 'genome') {
        title.textContent = M.label + '’s genome, all ' + nY + ' chromosomes';
        kicker.textContent = 'Figure 1 · Genome map · ' + M.label + '’s genome';
        mkick.textContent = 'The whole genome';
        mline.textContent = '';
      } else {
        var v0 = view(), band = bandName(ch, cam.c), gn = close_(lv) ? geneAt(ch, cam.c) : null;
        title.textContent = lv === 'chromosome' ? 'Chromosome ' + short(ch.n)
          : lv === 'region' ? 'Chromosome ' + short(ch.n) + (band ? ', band ' + band : '')
            : gn ? gn + ' and its neighbours' : 'Chromosome ' + short(ch.n) + ', ' + spanWords(v0[2]) + ' around ' + where(cam.c, v0[2]);
        kicker.textContent = 'Figure 1 · Genome map · ' + M.label;
        // Genome › chr13 › 13q13.1 › 32.30-32.42M: each step up is a button, the place itself is not.
        var sep = function () { return el('span', { 'class': 'sep', 'aria-hidden': 'true' }, ' › '); };
        if (lv === 'chromosome') { add(crumb, sep(), el('span', { 'aria-current': 'location' }, ch.n)); } else {
          var cb = el('button', { type: 'button' }, ch.n);
          cb.addEventListener('click', function () { openChrom(tgt.ci, undefined, undefined, true); });
          add(crumb, sep(), cb);
          var bi = bandIndex(ch, cam.c), bb = M.bands[ch.n];
          if (bi >= 0) {
            var bt = el('button', { type: 'button' }, band);
            bt.addEventListener('click', function () { openChrom(tgt.ci, (bb.s[bi] + bb.e[bi]) / 2, (bb.e[bi] - bb.s[bi]) * 1.1, true); });
            add(crumb, sep(), bt);
          }
          add(crumb, sep(), el('span', { 'aria-current': 'location' }, locusWords().replace(/^chr[^:]*:/, '')));
        }
        mkick.textContent = 'Under the line';
        mline.textContent = ch.n + ' · ' + where(cam.c, v0[2]) + (band ? ' · band ' + band : '');
      }
      var v = view();
      var inv = lv === 'genome' ? fs.length : fs.filter(function (r) { return r.chrom === ch.n && r.e > v[0] && r.s < v[1]; }).length;
      canvas.setAttribute('aria-label', lv === 'genome' ? 'The genome map: ' + nY + ' chromosomes with ' + fs.length + ' findings drawn'
        : 'Chromosome ' + short(ch.n) + ' open on the genome map, ' + where(v[0], v[2]) + ' to ' + where(v[1], v[2]) + ', ' +
          inv + (inv === 1 ? ' finding' : ' findings') + ' in view');
      narr.textContent = '';
      var shown = lv === 'genome' ? fs : fs.filter(function (r) { return r.chrom === ch.n && r.e > v[0] && r.s < v[1]; });
      var forms = (M.text.findings || {});
      if (lv === 'genome') {
        var form = fs.length === 0 ? forms.none : fs.length === 1 ? forms.one : forms.many;
        add(narr, runs((form || []).map(function (x) { return { t: x.t.replace('{n}', String(fs.length)), b: x.b }; })));
      } else if (M.text.chrom && M.text.chrom[ch.n] && rohLayer() === M.layers.roh) {
        // The chromosome's sentence is about its runs of homozygosity: a replay before their step ended has none.
        add(narr, runs(M.text.chrom[ch.n]));
      }
      var lay = M.layers.findings || {};
      if (!landed(lay)) {
        add(narr, el('p', {}, layerWords(lay)));
      } else if (shown.length && !(phone && lv === 'genome')) {
        // On a phone the genome level leaves the list to the sections under it: the page is long already.
        var list = shown.slice().sort(function (a, b) { return PRI[a.level] - PRI[b.level] || M.by[a.chrom] - M.by[b.chrom] || a.s - b.s; });
        var p = el('p', { 'class': 'gm-phrases' });
        var cap = phone ? 3 : 6;
        list.slice(0, cap).forEach(function (r, i) { if (i) { add(p, ' '); } add(p, r.phrase); });
        if (list.length > cap) { add(p, ' And ' + (list.length - cap) + ' more on the map.'); }
        add(narr, p);
      } else if (lv !== 'genome') {
        add(narr, el('p', {}, fs.some(function (r) { return r.chrom === ch.n; }) ? 'No finding in this stretch.' : noFindingHere()));
      }
      if (lv === 'genome') {
        var rl = rohLayer();
        var tail = rl === M.layers.roh ? runs(M.text.genome) : runs([{ t: layerWords(rl), b: false }]), lead = narr.children.length === 1 ? narr.firstChild : null;
        if (phone && lead) {
          // On a phone the count and the runs share one paragraph.
          add(lead, ' ');
          while (tail.firstChild) { lead.appendChild(tail.firstChild); }
        } else {
          add(narr, tail);
        }
        if (!phone) { add(narr, el('p', { 'class': 'dim' }, 'Click a chromosome to open it.')); }
      } else {
        var r = M.roh[ch.n], n = 0, longest = 0;
        if (r) { for (var i = first(r, v[0]); i < r.n && r.s[i] < v[1]; i++) { n++; longest = Math.max(longest, r.e[i] - r.s[i]); } }
        if (n && landed(rohLayer())) {
          add(narr, el('p', {}, 'Runs of homozygosity in view: ' + n + ', the longest ' + mb(longest) + ' Mb.'));
        }
        if (!phone) { add(narr, el('p', { 'class': 'dim' }, 'Drag to move along it, scroll or use + and − to zoom, / to find a gene, Esc to step back out.')); }
      }
      var nr = nearest();
      nearBtn.hidden = !nr;
      nearBtn.textContent = nr ? nr.words : '';
      drawCard();
      drawLayers();
      drawMarks();
      drawTable();
      stripWords();
    }
    function bandIndex(ch, bp) {
      var b = M.bands[ch.n];
      if (!b) { return -1; }
      for (var i = first(b, bp); i < b.n && b.s[i] <= bp; i++) { if (b.e[i] > bp) { return i; } }
      return -1;
    }
    function bandName(ch, bp) { var i = bandIndex(ch, bp); return i < 0 ? '' : short(ch.n) + M.bands[ch.n].o.name[i]; }
    function spanWords(span) {
      return span >= 1e6 ? mb(span) + ' Mb' : span >= 1e3 ? Math.round(span / 1e3) + ' kb' : Math.round(span) + ' bp';
    }
    function drawLayers() {
      layersBox.textContent = '';
      [['findings', 'Findings'], ['roh', 'Runs of homozygosity']].forEach(function (x) {
        var lay = x[0] === 'roh' ? rohLayer() : M.layers[x[0]];
        var row = el('div', { 'class': 'gm-layer', 'data-layer': x[0], 'data-state': lay ? lay.state : 'missing' });
        var dd = el('dd', {});
        if (!lay) { add(dd, 'Not on the map: this server’s map has no ' + x[1].toLowerCase() + ' layer.'); }
        else if (!landed(lay)) { add(dd, layerWords(lay)); }
        else {
          var n = x[0] === 'findings' ? drawn().length : lay.count;
          add(dd, n === 0 && lay.empty ? lay.empty : n + ' on the map. ');
          if (lay.note && n !== 0) { add(dd, lay.note); }
          if (lay.caveat) {
            var d = el('details', { 'class': 'gm-caveat' });
            add(d, el('summary', {}, 'Caveat'), el('p', {}, lay.caveat));
            add(dd, d);
          }
        }
        add(row, el('dt', {}, x[1] + (lay && lay.step ? ' · step ' + lay.step : '')), dd);
        add(layersBox, row);
      });
    }
    // The same marks keep their links, so the keyboard stays on the one it is on while the map moves.
    var marksKey = '';
    function drawMarks() {
      var lv = level(), ch = chrom(), v = view();
      var list = drawn().filter(function (r) { return lv === 'genome' || (r.chrom === ch.n && r.e > v[0] && r.s < v[1]); })
        .sort(function (a, b) { return M.by[a.chrom] - M.by[b.chrom] || a.s - b.s; });
      var k = list.map(function (r) { return r.id + '\u0000' + r.title + '\u0000' + r.href; }).join('\u0001');
      if (k === marksKey && marks.childNodes.length === list.length) { return; }
      marksKey = k;
      marks.textContent = '';
      list.forEach(function (r) {
          var a = el('a', { href: r.href, 'data-id': r.id, 'data-level': r.level, 'data-chrom': r.chrom,
            'data-s': String(r.s), 'data-e': String(r.e) }, (U.LEVEL_WORDS[r.level] || r.level) + ': ' + r.title + ', ' + r.chrom);
          // Tab through the marks: the one with the keyboard gets a ring on the map and its card.
          a.addEventListener('focus', function () { focusId = r.id; pick(r.id); });
          a.addEventListener('blur', function () { if (focusId === r.id) { focusId = null; cacheKey = ''; schedule(); } });
          add(marks, add(el('li', {}), a));
        });
    }
    function drawKey() {
      keyBox.textContent = '';
      add(keyBox, el('div', { 'class': 'fig-label' }, 'Key'));
      ['act', 'medicine', 'review', 'carrier', 'info'].forEach(function (lv) { add(keyBox, U.level(lv)); });
      add(keyBox, add(el('span', { 'class': 'gm-k-roh' }), el('span', { 'class': 'sw', 'aria-hidden': 'true' }), 'Run of homozygosity, 5 Mb or more'));
      add(keyBox, add(el('span', { 'class': 'gm-k-mask' }), el('span', { 'class': 'sw', 'aria-hidden': 'true' }), 'Hatched: centromere or heterochromatin, few reads map there'));
    }

    // --- the card, Nearest and the search: each reads drawn(), so a closed result is never offered ---
    function rowOf(id) { return id ? drawn().filter(function (r) { return r.id === id; })[0] || null : null; }
    // The card of a finding on the map (or none). It lives in the page's memory only, never in the address.
    function pick(id) {
      sel = rowOf(id) ? id : null;
      drawCard();
      cacheKey = '';
      schedule();
    }
    function drawCard() {
      var r = rowOf(sel);
      card.textContent = '';
      if (!r) { sel = null; card.hidden = true; return; }
      card.hidden = false;
      var where_ = r.chrom + ':' + (r.pos ? comma(r.pos) : comma(r.s + 1) + '-' + comma(r.e));
      add(card, el('div', { 'class': 'fig-label' }, U.LEVEL_WORDS[r.level] || r.level));
      add(card, el('h3', { 'class': 'gm-card-title' }, r.title));
      var facts = el('p', { 'class': 'gm-card-facts mono' }, (r.gene ? r.gene + ' · ' : '') + where_ + (r.hgvs ? ' · ' + r.hgvs : ''));
      add(card, facts);
      if (r.phrase) { add(card, el('p', { 'class': 'gm-card-phrase' }, r.phrase)); }
      var foot = el('p', { 'class': 'gm-card-foot' });
      add(foot, el('a', { href: r.href }, 'The finding’s page ›'), ' ', button('Close the card', function () { pick(null); canvas.focus({ preventScroll: true }); }));
      add(card, foot);
    }
    // The finding nearest the view when none of this chromosome's is in it: the closest one on this
    // chromosome, else the first on the next chromosome that has one, in reading order. Null when one is in view.
    function nearest() {
      if (!M || level() === 'genome') { return null; }
      var ch = chrom(), v = view(), fs = drawn(), best = null, d = Infinity, side = 0;
      for (var i = 0; i < fs.length; i++) {
        var r = fs[i];
        if (r.chrom !== ch.n) { continue; }
        if (r.e > v[0] && r.s < v[1]) { return null; }
        var dr = r.s >= v[1] ? r.s - v[1] : v[0] - r.e;
        if (dr < d) { d = dr; best = r; side = r.s >= v[1] ? 1 : -1; }
      }
      var name;
      if (best) {
        name = (best.gene || best.title) + ' at ' + mb((best.s + best.e) / 2) + ' Mb';
        return { r: best, side: side, words: side < 0 ? '‹ Nearest: ' + name : 'Nearest: ' + name + ' ›' };
      }
      for (var k = 1; k < M.chroms.length && !best; k++) {
        var n = M.chroms[(cam.ci + k) % M.chroms.length].n;
        fs.forEach(function (x) { if (x.chrom === n && (!best || x.s < best.s)) { best = x; } });
      }
      if (!best) { return null; }
      return { r: best, side: 1, words: 'Nearest: ' + (best.gene || best.title) + ' on ' + best.chrom + ' ›' };
    }
    // Fly to a finding at the gene level and open its card.
    function goTo(r) {
      var span = Math.max((r.e - r.s) * 1.5, GENE_STOP * openW());
      openChrom(M.by[r.chrom], (r.s + r.e) / 2, span, true);
      pick(r.id);
    }
    function goNearest() { var n = nearest(); if (n) { goTo(n.r); } }

    // Search: the findings drawn now (gene, id, title), the gene names (the names file, fetched on the
    // first key), the band names, and a place written as the address writes it.
    var found = [], active = -1;
    function loadNames() {
      if (names || !M.gencode) { return; }
      names = { state: 'loading' };
      fetchGz(genesUrl('names')).then(function (d) {
        names = { state: 'ready', name: d.name || [], up: (d.name || []).map(function (x) { return String(x).toUpperCase(); }),
          chrom: d.chrom || [], s: d.s || [], e: d.e || [] };
        if (document.activeElement === findBox) { runFind(); }
      }, function (err) {
        names = { state: 'failed', note: String((err && err.message) || err).replace(/\.$/, '') };
        if (document.activeElement === findBox) { runFind(); }
      });
    }
    function findAll(q) {
      var out = [], Q = q.toUpperCase(), bare = Q.replace(/^CHR/, '');
      if (!q) { return out; }
      var p = /^chr/i.test(q) && /[:]/.test(q) ? parsePlace(q.replace(/\s+/g, '')) : null;
      if (p && p.k) { out.push({ kind: 'place', words: 'Go to ' + q, p: p }); }
      var fs = drawn();
      fs.filter(function (r) {
        return (r.gene && r.gene.toUpperCase().indexOf(Q) === 0) || r.id.toUpperCase().indexOf(Q) === 0 ||
          (Q.length >= 3 && String(r.title).toUpperCase().indexOf(Q) >= 0);
      }).sort(function (a, b) { return PRI[a.level] - PRI[b.level] || M.by[a.chrom] - M.by[b.chrom] || a.s - b.s; }).slice(0, 3).forEach(function (r) {
        out.push({ kind: 'finding', words: (r.gene ? r.gene + ' · ' : '') + (U.LEVEL_WORDS[r.level] || r.level) + ': ' + r.title, r: r });
      });
      if (names && names.state === 'ready') {
        var exact = [], pre = [];
        for (var i = 0; i < names.up.length && pre.length < 5; i++) {
          var u = names.up[i];
          if (u.indexOf(Q) !== 0 || !(names.chrom[i] in M.by)) { continue; }
          (u === Q ? exact : pre).push(i);
        }
        exact.concat(pre).slice(0, 5).forEach(function (j) {
          out.push({ kind: 'gene', words: names.name[j] + ' · gene · ' + names.chrom[j], c: names.chrom[j], s: names.s[j], e: names.e[j] });
        });
      }
      var nb = 0;
      M.chroms.forEach(function (ch) {
        var b = M.bands[ch.n];
        for (var i2 = 0; b && i2 < b.n && nb < 2; i2++) {
          var full = short(ch.n) + b.o.name[i2];
          if (bare.length >= 2 && full.toUpperCase().indexOf(bare) === 0) {
            nb++;
            out.push({ kind: 'band', words: full + ' · band', c: ch.n, s: b.s[i2], e: b.e[i2] });
          }
        }
      });
      return out.slice(0, 8);
    }
    function runFind() {
      var q = findBox.value.trim();
      if (q) { loadNames(); }
      found = M ? findAll(q) : [];
      active = found.length ? 0 : -1;
      drawFound(q);
    }
    function drawFound(q) {
      findList.textContent = '';
      var note = '';
      if (q && names && names.state === 'loading') { note = 'Looking through the gene names…'; }
      else if (q && names && names.state === 'failed') { note = 'Gene names: ' + names.note + '.'; }
      else if (q && !M.gencode) { note = 'Gene names need step 08’s GENCODE file.'; }
      found.forEach(function (x, i) {
        var li = el('li', { role: 'option', id: 'gm-find-o' + i, 'data-kind': x.kind, 'aria-selected': String(i === active) }, x.words);
        li.addEventListener('mousedown', function (e) { e.preventDefault(); choose(i); });
        add(findList, li);
      });
      if (q && !found.length && !note) { note = 'Nothing on the map by that name.'; }
      if (note) { add(findList, el('li', { 'class': 'gm-find-note', role: 'presentation' }, note)); }
      findList.hidden = !q;
      findBox.setAttribute('aria-expanded', String(!!q && found.length > 0));
      if (active >= 0) { findBox.setAttribute('aria-activedescendant', 'gm-find-o' + active); } else { findBox.removeAttribute('aria-activedescendant'); }
    }
    function closeFind() { found = []; active = -1; findList.textContent = ''; findList.hidden = true; findBox.setAttribute('aria-expanded', 'false'); findBox.removeAttribute('aria-activedescendant'); }
    function choose(i) {
      var x = found[i];
      if (!x || !M) { return; }
      closeFind();
      if (x.kind === 'finding') { goTo(x.r); }
      else if (x.kind === 'place') { openChrom(x.p.ci, x.p.s === undefined ? undefined : (x.p.s + x.p.e) / 2, x.p.s === undefined ? undefined : x.p.e - x.p.s, true); }
      else { openChrom(M.by[x.c], (x.s + x.e) / 2, Math.max((x.e - x.s) * (x.kind === 'gene' ? 1.3 : 1.1), 60000), true); }
      canvas.focus({ preventScroll: true });
    }
    findBox.addEventListener('input', runFind);
    findBox.addEventListener('focus', function () { if (findBox.value.trim()) { runFind(); } });
    findBox.addEventListener('blur', function () { window.setTimeout(closeFind, 150); });
    findBox.addEventListener('keydown', function (e) {
      if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
        if (found.length) { active = (active + (e.key === 'ArrowDown' ? 1 : found.length - 1)) % found.length; drawFound(findBox.value.trim()); }
        e.preventDefault();
      } else if (e.key === 'Enter') {
        if (found.length) { choose(active < 0 ? 0 : active); }
        e.preventDefault();
      } else if (e.key === 'Escape') {
        findBox.value = '';
        closeFind();
        canvas.focus({ preventScroll: true });
        e.preventDefault();
      }
    });

    // Paint the copies: none on a person's page. Origin is a family run's answer, so it says why it is off.
    function inFamily() {
      var st = window.PGP && window.PGP.state;
      if (!st || !st.plugins) { return null; }
      var fam = st.plugins.family;
      if (fam && fam.unknown) { return null; }  // a family report the server could not read: it cannot say
      return !!(fam && fam.kind === 'member' && (fam.families || []).length);
    }
    function drawPaint() {
      paintBox.textContent = '';
      var fam = inFamily();
      var why = fam ? 'Origin: which parent each stretch came from is painted on the family page.'
        : fam === false ? 'Origin needs a family run with a sequenced parent: this run is in none.'
          : 'Origin needs a family run with a sequenced parent.';
      add(paintBox, el('div', { 'class': 'fig-label' }, 'Paint the copies'));
      var seg = el('div', { 'class': 'gm-seg' });
      add(seg, el('button', { type: 'button', 'data-paint': 'none', 'aria-pressed': 'true' }, 'None'));
      add(seg, el('button', { type: 'button', 'data-paint': 'origin', 'aria-pressed': 'false', 'aria-disabled': 'true',
        'aria-describedby': 'gm-paint-why' }, 'Origin'));
      add(paintBox, seg, el('p', { 'class': 'gm-paint-why', id: 'gm-paint-why' }, why));
    }

    // --- the closed results: one dashed chip per kind, and the opened ones ---
    function drawChips() {
      chips.textContent = '';
      if (!M) { return; }
      var acts = closedActs(), gates = [], opened = [];
      M.closable.forEach(function (id) {
        // A gated result (APOE) has a chip only once the page shows its card (a replay before its step has none).
        if (M.gated[id] && opts.now && !opts.now(id)) { return; }
        var s = stateOf(id);
        if (s === 'closed') { if (M.gated[id]) { gates.push(id); } } else if (s === 'opened') { opened.push(id); }
      });
      // Why it is closed: said beside the count, left out on a phone (the card says it too).
      function why() { return el('span', { 'class': 'gm-why' }, ' · Not drawn until opened on purpose'); }
      if (acts.length) {
        var c = el('div', { 'class': 'gm-chip closed', 'data-closed': 'act' });
        add(c, add(el('span', {}, acts.length === 1 ? '1 result closed' : acts.length + ' results closed'), why()));
        add(c, confirming === 'act' ? confirmBox(acts[0], false) : button('Open', function () { confirming = 'act'; drawChips(); keyTo('.gm-yes'); }));
        add(chips, c);
      }
      gates.forEach(function (id) {
        var g = M.gated[id], name = g.gene || 'This result';
        var c = el('div', { 'class': 'gm-chip closed', 'data-closed': 'gate' });
        add(c, add(el('span', {}, name + ' closed'), why()));
        add(c, confirming === id ? confirmBox(id, true) : button('Open ' + name, function () { confirming = id; drawChips(); keyTo('.gm-yes'); }));
        add(chips, c);
      });
      opened.forEach(function (id) {
        var r = M.findings.filter(function (x) { return x.id === id; })[0];
        var c = el('div', { 'class': 'gm-chip opened', 'data-opened': '1' });
        add(c, el('span', {}, 'Opened: ' + (r ? r.title : 'a result to act on')));
        add(c, button('Close it again', function () { confirming = null; opts.close(id); keyTo(); }));
        add(chips, c);
      });
      // A result the chip opened that has no place: named only once it is open, with a link to its card.
      inSection = inSection.filter(function (id) { return stateOf(id) !== 'closed' && opts.title && opts.title(id); });
      inSection.forEach(function (id) {
        var c = el('div', { 'class': 'gm-chip opened', 'data-opened': 'section' });
        add(c, add(el('span', {}, 'Opened in its section: '), el('a', { href: '#f/' + id }, opts.title(id))));
        add(c, button('Close it again', function () { confirming = null; opts.close(id); keyTo(); }));
        add(chips, c);
      });
    }
    // Where the keyboard goes once a chip's button rebuilt the chips: the chips' first button, else the map.
    function keyTo(sel) {
      var b = chips.querySelector(sel || 'button');
      (b || canvas).focus({ preventScroll: true });
    }
    function confirmBox(id, gated) {
      var box = el('span', { 'class': 'gm-confirm' });
      add(box, el('span', {}, gated ? 'It is kept behind its warning on every page. Open it here and in its section?'
        : 'A result to act on in adult life is usually not looked at before the child can choose. Open it here and in its section?'));
      add(box, button('Show the result', function () {
        confirming = null;
        if (!gated && M.closable.indexOf(id) < 0 && inSection.indexOf(id) < 0) { inSection.push(id); }
        opts.open(id);  // the page draws again, and the map's chips with it (refresh)
        keyTo();
      }, 'gm-yes'));
      add(box, button('Keep it closed', function () { confirming = null; drawChips(); keyTo(); }, 'gm-no'));
      return box;
    }
    function button(text, fn, cls) {
      var b = el('button', { type: 'button', 'class': cls || '' }, text);
      b.addEventListener('click', fn);
      return b;
    }

    // --- input ---
    function at(e) { var r = canvas.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; }
    function hitAt(x, y) {
      var best = null;
      for (var i = hits.length - 1; i >= 0; i--) {
        var h = hits[i];
        if (h.kind === 'finding') { if (Math.abs(h.x - x) <= h.rad && Math.abs(h.y - y) <= h.rad) { return h; } }
        else if (x >= h.x0 && x <= h.x1 && y >= h.y0 && y <= h.y1) {
          if (!best || (h.kind === 'masked' && best.kind === 'band') || (h.kind === 'roh')) { best = h; }
        }
      }
      return best;
    }
    function tipFor(h) {
      var ch = M.chroms[cam.ci];
      if (h.kind === 'finding') { return [h.r.phrase, 'Click for its card']; }
      if (h.kind === 'gene') {
        return [h.name + (h.type ? ' · ' + h.type.replace(/_/g, ' ') : ''),
          ch.n + ':' + comma(h.s + 1) + '-' + comma(h.e) + (h.strand === '+' || h.strand === '-' ? ' · ' + h.strand + ' strand' : '') + ' · click to zoom to it'];
      }
      if (h.kind === 'nearest') { return null; }
      if (h.kind === 'band') { return ['Band ' + short(ch.n) + (h.name || ''), ch.n + ':' + comma(h.s + 1) + '-' + comma(h.e)]; }
      if (h.kind === 'masked') {
        return ['Not readable here · ' + ch.n + ':' + comma(h.s + 1) + '-' + comma(h.e),
          'Centromere, heterochromatin or acrocentric arm: few reads map, so layers are cut here. Not missing data.'];
      }
      if (h.kind === 'roh') {
        return [mb(h.e - h.s) + ' Mb run of homozygosity' + (h.art ? ' (' + h.art + ': not a real run)' : ''), ch.n + ':' + comma(h.s + 1) + '-' + comma(h.e)];
      }
      if (h.kind === 'row') { return null; }
      return null;
    }
    // A hidden tooltip keeps no words: it may have named a result that is closed now.
    function hideTip() { tip.hidden = true; tip.textContent = ''; }
    function showTip(h, x, y) {
      var t = h && tipFor(h);
      if (!t) { hideTip(); return; }
      tip.textContent = '';
      add(tip, el('b', {}, t[0]));
      if (t[1]) { add(tip, el('span', {}, t[1])); }
      tip.hidden = false;
      tip.style.setProperty('left', Math.min(x + 12, Math.max(0, CW - 260)) + 'px');
      tip.style.setProperty('top', (y + 14) + 'px');
    }
    var drag = null;
    canvas.addEventListener('pointermove', function (e) {
      if (!M) { return; }
      var p = at(e);
      if (drag && level() !== 'genome') {
        var v = view();
        tgt.c = drag.c - (p[0] - drag.x) / plotW() * v[2];
        clampTgt();
        cam.c = tgt.c;
        drag.moved = drag.moved || Math.abs(p[0] - drag.x) > 3;
        changed(false);
        return;
      }
      var h = hitAt(p[0], p[1]);
      var row = level() === 'genome' && h ? (h.kind === 'row' ? h.i : h.kind === 'finding' ? M.by[h.r.chrom] : -1) : -1;
      if (row !== hover.row) { hover.row = row; schedule(); }
      hover.hit = h;
      canvas.style.setProperty('cursor', h && /^(finding|row|gene|nearest)$/.test(h.kind) ? 'pointer' : level() === 'genome' ? 'default' : 'grab');
      showTip(h, p[0], p[1]);
    });
    canvas.addEventListener('pointerleave', function () { hover.row = -1; hideTip(); schedule(); });
    canvas.addEventListener('pointerdown', function (e) {
      if (!M || level() === 'genome' || e.button !== 0) { return; }
      drag = { x: at(e)[0], c: tgt.c, moved: false };
      if (canvas.setPointerCapture) { canvas.setPointerCapture(e.pointerId); }
    });
    canvas.addEventListener('pointerup', function (e) {
      var was = drag;
      drag = null;
      if (!M || (was && was.moved)) { if (was && was.moved) { writeHash(false); } return; }
      var p = at(e), h = hitAt(p[0], p[1]);
      if (h && h.kind === 'finding') { pick(h.r.id); return; }
      if (h && h.kind === 'nearest') { goNearest(); return; }
      if (h && h.kind === 'gene') { openChrom(tgt.ci, (h.s + h.e) / 2, Math.max((h.e - h.s) * 1.3, 60000), true); return; }
      if (level() === 'genome' && h && h.kind === 'row') { openChrom(h.i, undefined, undefined, true); }
    });
    canvas.addEventListener('wheel', function (e) {
      if (!M || level() === 'genome') { return; }
      e.preventDefault();
      var p = at(e), v = view(), ch = chrom();
      var bp = v[0] + (p[0] - gut()) / plotW() * v[2];
      var f = Math.exp(clamp(e.deltaY, -200, 200) * 0.0025);
      if (f > 1 && Math.exp(tgt.ls) >= ch.len * 0.999) { closeChrom(true); return; }
      var span = clamp(Math.exp(tgt.ls) * f, minSpan(ch.len), ch.len);
      tgt.c = bp - (bp - tgt.c) * span / Math.exp(tgt.ls);
      tgt.ls = Math.log(span);
      clampTgt();
      changed(false);
    }, { passive: false });
    // The map's keys: arrows pan a tenth of the view (Shift: a whole view), + and − zoom 2× at the centre,
    // [ and ] step through the chromosomes, 0 is the whole genome, and Esc closes the card, then steps out one
    // level at a time. Tab is left alone: it steps through the marks' links, then out of the map.
    canvas.addEventListener('keydown', function (e) {
      if (!M || e.altKey || e.ctrlKey || e.metaKey) { return; }
      var open = level() !== 'genome', span = Math.exp(tgt.ls);
      var done = true;
      if (e.key === 'Escape') {
        if (sel) { pick(null); } else if (open) { stepOut(); } else { done = false; }
      } else if (e.key === '0') { if (open) { closeChrom(true); } else { done = false; } }
      else if (e.key === '[' || e.key === ']') {
        var i = (open ? tgt.ci : (e.key === ']' ? -1 : M.chroms.length)) + (e.key === ']' ? 1 : -1);
        openChrom((i + M.chroms.length) % M.chroms.length, undefined, undefined, true);
      } else if (open && (e.key === '+' || e.key === '=' || e.key === '-' || e.key === '−')) {
        zoomBy(e.key === '-' || e.key === '−' ? 2 : 0.5);
      } else if (open && (e.key === 'ArrowLeft' || e.key === 'ArrowRight')) {
        tgt.c += (e.key === 'ArrowLeft' ? -1 : 1) * (e.shiftKey ? 1 : 0.1) * span;
        clampTgt();
        changed(false);
      } else if (!open && e.key === 'Enter') { openChrom(tgt.ci, undefined, undefined, true); }
      else { done = false; }
      if (done) { e.preventDefault(); }
    });
    // Esc, once no card is open: reads to the gene, gene to the region, region to the chromosome, chromosome to the genome.
    function stepOut() {
      var lv = level();
      if (lv === 'chromosome') { closeChrom(true); }
      else if (lv === 'region') { openChrom(tgt.ci, undefined, undefined, true); }
      else { openChrom(tgt.ci, tgt.c, (lv === 'reads' ? GENE_STOP : REGION_STOP) * openW(), true); }
    }
    // The keys of the whole page while the map is on it: / finds, L opens the layers. Never inside a text box,
    // never with Cmd, Ctrl or Alt.
    document.addEventListener('keydown', function (e) {
      if (!M || failure || e.altKey || e.ctrlKey || e.metaKey || !fig.isConnected) { return; }
      var t = e.target, tag = t && t.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT' || (t && t.isContentEditable)) { return; }
      if (e.key === '/') { e.preventDefault(); findBox.focus(); }
      else if (e.key === 't' || e.key === 'T') { e.preventDefault(); inview.open = !inview.open; }
      else if (e.key === 'l' || e.key === 'L') {
        e.preventDefault();
        more.open = true;
        var b = paintBox.querySelector('button');
        if (b) { b.focus(); }
      }
    });
    function zoomBy(f) {
      if (!M || level() === 'genome') { return; }
      tgt.ls = Math.log(Math.exp(tgt.ls) * f);
      clampTgt();
      changed(false);
    }
    zoomOut.addEventListener('click', function () { zoomBy(2); });
    zoomIn.addEventListener('click', function () { zoomBy(0.5); });
    stopG.addEventListener('click', function () { closeChrom(true); });
    // The other stops keep the place in view; from the genome they go to the first finding by level (act
    // first), else to the middle of chromosome 1.
    function stopTo(lv) {
      var ci = tgt.ci, c = tgt.c;
      if (level() === 'genome') {
        var f = drawn().slice().sort(function (a, b) { return PRI[a.level] - PRI[b.level] || M.by[a.chrom] - M.by[b.chrom] || a.s - b.s; })[0];
        ci = f ? M.by[f.chrom] : 0;
        c = f ? (f.s + f.e) / 2 : M.chroms[ci].len / 2;
      }
      if (lv === 'chromosome') { openChrom(ci, undefined, undefined, true); return; }
      openChrom(ci, c, (lv === 'region' ? REGION_STOP : GENE_STOP) * openW(), true);
    }
    stopC.addEventListener('click', function () { stopTo('chromosome'); });
    stopR.addEventListener('click', function () { stopTo('region'); });
    stopN.addEventListener('click', function () { stopTo('gene'); });
    nearBtn.addEventListener('click', function () { goNearest(); });
    window.addEventListener('hashchange', function () { readHash(false); });
    window.addEventListener('popstate', function () { readHash(false); });
    if (window.ResizeObserver) { new window.ResizeObserver(function () { if (M) { size(); clampTgt(); schedule(); } }).observe(stage); }

    // --- the whole-genome strip: where the open view sits in the genome ---
    // Every chromosome end to end on one axis, 5 Mb apart, a tick for each finding drawn (a closed one is never
    // drawn, so never ticked), and the view as a brass brush at least 3 px wide with the rest dimmed. It shows
    // while a chromosome is open. On a desktop a click opens a chromosome whole, a drag on the brush pans and a
    // drag elsewhere brushes a new view on the chromosome where it started, clamped to its ends; on a phone the
    // strip is a picture only (chromosome 21 is a few pixels wide there).
    var STRIP_GAP = 5e6, SW = 0, SH = 0, brush = null;
    function stripGeo() {
      var x0 = phone ? 64 : 128, x1 = Math.max(x0 + 50, CW - pad()), off = [], at = 0;
      M.chroms.forEach(function (c, i) { off[i] = at; at += c.len + STRIP_GAP; });
      return { x0: x0, x1: x1, off: off, per: (x1 - x0) / (at - STRIP_GAP) };
    }
    function sx(G, i, bp) { return G.x0 + (G.off[i] + bp) * G.per; }
    // The chromosome under x on the strip, -1 in a gap.
    function stripChrom(G, x) {
      for (var i = 0; i < M.chroms.length; i++) { if (x >= sx(G, i, 0) - 1 && x <= sx(G, i, M.chroms[i].len) + 1) { return i; } }
      return -1;
    }
    // The brush: the view, or the span being brushed.
    function stripBrush(G) {
      var ch = chrom(), v = view(), a, z;
      if (brush && brush.mode === 'brush' && brush.moved) {
        a = sx(G, brush.ci, Math.min(brush.a, brush.z));
        z = sx(G, brush.ci, Math.max(brush.a, brush.z));
      } else {
        a = sx(G, cam.ci, clamp(v[0], 0, ch.len));
        z = sx(G, cam.ci, clamp(v[1], 0, ch.len));
      }
      if (z - a < 3) { var m = (a + z) / 2; a = m - 1.5; z = m + 1.5; }
      return [a, z];
    }
    function drawStrip() {
      var on = Math.max(cam.k, tgt.k) >= 0.5;
      if (strip.hidden === on) { strip.hidden = !on; }
      if (!on) { strip.pgpDrawn = []; return; }
      var h = phone ? 12 : 34;
      if (CW !== SW || h !== SH || strip.width !== Math.round(CW * DPR)) {
        SW = CW; SH = h;
        strip.width = Math.round(CW * DPR);
        strip.height = Math.round(h * DPR);
        strip.style.setProperty('width', CW + 'px');
        strip.style.setProperty('height', h + 'px');
      }
      var c = strip.getContext('2d'), G = stripGeo(), yc = phone ? 6 : 15, bh = phone ? 5 : 8, th = phone ? 10 : 12;
      c.setTransform(1, 0, 0, 1, 0, 0);
      c.clearRect(0, 0, strip.width, strip.height);
      c.setTransform(DPR, 0, 0, DPR, 0, 0);
      c.globalAlpha = clamp((cam.k - 0.5) * 2, 0, 1);
      c.fillStyle = INK3;
      c.font = (phone ? '600 8px ' : '600 10.5px ') + SANS;
      c.textAlign = 'right';
      c.textBaseline = 'middle';
      c.fillText(phone ? 'GENOME' : 'WHOLE GENOME', G.x0 - (phone ? 6 : 12), yc);
      var nf = '9.5px ' + SANS;
      M.chroms.forEach(function (ch, i) {
        var a = sx(G, i, 0), z = sx(G, i, ch.len);
        c.fillStyle = '#2a2c31';
        c.beginPath(); rounded(c, a, yc - bh / 2, Math.max(1, z - a), bh, bh / 2); c.fill();
        var nm = short(ch.n);
        if (!phone && z - a > measure(nf, nm) + 2) {
          c.fillStyle = INK3;
          c.font = nf;
          c.textAlign = 'center';
          c.textBaseline = 'top';
          c.fillText(nm, (a + z) / 2, yc + 10);
        }
      });
      var ticked = [];
      drawn().forEach(function (r) {
        if (!drawsAt('findings', r, 'genome')) { return; }
        var x = sx(G, M.by[r.chrom], (r.s + r.e) / 2);
        c.fillStyle = LV[r.level] || LV.info;
        c.fillRect(x - 1, yc - th / 2, 2, th);
        ticked.push(mark('findings', r.chrom, r.s, r.e, r.id));
      });
      var B = stripBrush(G);
      c.fillStyle = 'rgba(14,15,16,.55)';
      c.fillRect(0, 0, Math.max(0, B[0]), h);
      c.fillRect(B[1], 0, Math.max(0, CW - B[1]), h);
      c.beginPath(); rounded(c, B[0], yc - th / 2 - 2, B[1] - B[0], th + 4, 3);
      c.fillStyle = 'rgba(200,171,120,.10)';
      c.fill();
      c.strokeStyle = BRASS;
      c.lineWidth = 1.5;
      c.stroke();
      if (brush && brush.mode === 'brush' && brush.moved) {
        // While brushing, the size of the span in brass beside the brush.
        var words = spanWords(Math.abs(brush.z - brush.a)), font = '600 10.5px ' + MONO, w = measure(font, words);
        c.textAlign = 'left';
        c.textBaseline = 'middle';
        label(c, words, B[1] + 6 + w < CW ? B[1] + 6 : B[0] - 6 - w, yc, font, BRASS);
      }
      c.globalAlpha = 1;
      strip.pgpDrawn = ticked;
    }
    function stripX(e) { return e.clientX - strip.getBoundingClientRect().left; }
    strip.addEventListener('pointerdown', function (e) {
      if (!M || phone || e.button !== 0 || level() === 'genome') { return; }
      var G = stripGeo(), x = stripX(e), B = stripBrush(G), inside = x >= B[0] - 3 && x <= B[1] + 3;
      var i = inside ? cam.ci : stripChrom(G, x);
      brush = { mode: inside ? 'pan' : 'brush', x: x, c: tgt.c, ci: i, moved: false, a: 0, z: 0 };
      if (i >= 0) { brush.a = brush.z = clamp((x - sx(G, i, 0)) / G.per, 0, M.chroms[i].len); }
      if (strip.setPointerCapture) { strip.setPointerCapture(e.pointerId); }
      e.preventDefault();
    });
    strip.addEventListener('pointermove', function (e) {
      if (!M || phone || level() === 'genome') { return; }
      var G = stripGeo(), x = stripX(e);
      if (!brush) {
        var B = stripBrush(G);
        strip.style.setProperty('cursor', x >= B[0] - 3 && x <= B[1] + 3 ? 'grab' : stripChrom(G, x) >= 0 ? 'pointer' : 'default');
        return;
      }
      brush.moved = brush.moved || Math.abs(x - brush.x) > 3;
      if (!brush.moved) { return; }
      if (brush.mode === 'pan') {
        tgt.c = brush.c + (x - brush.x) / G.per;
        clampTgt();
        cam.c = tgt.c;
        changed(false);
      } else if (brush.ci >= 0) {
        brush.z = clamp((x - sx(G, brush.ci, 0)) / G.per, 0, M.chroms[brush.ci].len);
        schedule();
      }
    });
    strip.addEventListener('pointerup', function (e) {
      var b0 = brush;
      brush = null;
      if (!M || phone || !b0) { return; }
      if (!b0.moved) {
        var i = stripChrom(stripGeo(), stripX(e));
        if (i >= 0) { openChrom(i, undefined, undefined, true); }
      } else if (b0.mode === 'pan') {
        writeHash(false);
      } else if (b0.ci >= 0) {
        var a = Math.min(b0.a, b0.z), z = Math.max(b0.a, b0.z);
        openChrom(b0.ci, (a + z) / 2, Math.max(z - a, MIN_SPAN), true);
      }
      schedule();
    });
    strip.addEventListener('pointercancel', function () { brush = null; schedule(); });
    // The strip's words for a screen reader, on settle: where the view is.
    function stripWords() {
      var lv = level(), v = view(), ch = chrom();
      strip.setAttribute('aria-label', lv === 'genome' ? 'The whole genome' : 'Where the view is in the whole genome: chromosome ' +
        short(ch.n) + ', ' + where(Math.max(0, v[0]), v[2]) + ' to ' + where(Math.min(ch.len, v[1]), v[2]));
    }

    // --- the In view table: every mark the map draws one by one in view ---
    // A real table under the map, the screen reader's path and the overflow for the names the canvas has no room
    // for. Its rows are what drawsAt() draws one by one in view, from what the page shows now (drawn(): a closed
    // result is never a row); a lane drawn as a density adds one line instead. Sorted by place with the findings
    // first; a header sorts by its column, again the other way. 12 rows, then Show all (a scrolled window past 200
    // rows, at most 5,000). It follows the map on settle, and its count is the one thing it announces. A row's
    // button flies the map to the feature, with a finding's card.
    var TCOLS = [['kind', 'Kind'], ['what', 'What'], ['where', 'Where'], ['size', 'Size'], ['who', 'Who'], ['copy', 'Copy'],
      ['source', 'Source']];
    var T_SHORT = 12, T_FULL = 200, T_MAX = 5000, ROW_PX = 30;
    var tSort = { col: 'where', dir: 1 }, showAll = false, announced = -1, tRows = [], selKey = null, tSig = '';
    var inCount = el('span', { 'class': 'gm-inview-n' });
    var inSum = add(el('summary', {}), el('span', { 'class': 'gm-inview-t' }, 'In view'), ' · ', inCount);
    var inLive = el('p', { 'class': 'gm-inview-live', 'aria-live': 'polite' });
    var inWrap = el('div', { 'class': 'gm-inview-wrap' });
    var inDense = el('div', { 'class': 'gm-inview-dense' });
    var inFoot = el('div', { 'class': 'gm-inview-foot' });
    add(inview, inSum, inLive, inWrap, inDense, inFoot);
    // Open beside a desktop map, folded on a phone (the page is long already). T opens and folds it.
    inview.open = !(window.matchMedia && window.matchMedia('(max-width: 700px)').matches);
    var COPY = { het: 'one copy (unphased)', hom: 'both copies', hemi: 'the only copy' };
    function sizeWords(bp) {
      if (bp < 1000) { return comma(bp) + ' bp'; }
      if (bp < 1e6) { return (bp / 1e3).toFixed(bp < 1e4 ? 1 : 0) + ' kb'; }
      return mb(bp) + ' Mb';
    }
    function placeWords(c, s, e, pos) {
      if (pos) { return c + ':' + comma(pos); }
      return c + ':' + (e - s === 1 ? comma(e) : comma(s + 1) + '-' + comma(e));
    }
    function stepWords(steps) {
      steps = steps || [];
      return steps.length ? (steps.length === 1 ? 'Step ' : 'Steps ') + steps.join(', ') : '';
    }
    function findingRow(r) {
      return { key: mark('findings', r.chrom, r.s, r.e, r.id), layer: 'findings', r: r, ci: M.by[r.chrom], s: r.s, e: r.e,
        rank: PRI[r.level], kind: U.LEVEL_WORDS[r.level] || r.level, what: r.title, where: placeWords(r.chrom, r.s, r.e, r.pos),
        size: r.e - r.s, who: M.label, copy: COPY[r.zyg] || 'not known', source: stepWords(r.steps) };
    }
    function rohRow(ch, R, i, lay) {
      var art = M.art[R.o.art[i]] || '';
      return { key: mark('roh', ch.n, R.s[i], R.e[i], ''), layer: 'roh', ci: M.by[ch.n], s: R.s[i], e: R.e[i], rank: 10,
        kind: 'Run of homozygosity', what: art ? 'Not a real run: ' + art : R.o.long[i] ? 'A long run, 5 Mb or more' : 'A run under 5 Mb',
        where: placeWords(ch.n, R.s[i], R.e[i]), size: R.e[i] - R.s[i], who: M.label, copy: art ? '' : 'both copies alike',
        source: lay.step ? 'Step ' + lay.step : '' };
    }
    function geneRow(ch, G, j) {
      var type = G.o.type[j] || '';
      return { key: mark('genes', ch.n, G.s[j], G.e[j], G.o.name[j]), layer: 'genes', ci: M.by[ch.n], s: G.s[j], e: G.e[j], rank: 20,
        kind: 'Gene', what: G.o.name[j] + (type && type !== 'protein_coding' ? ' · ' + type.replace(/_/g, ' ') : ''),
        where: placeWords(ch.n, G.s[j], G.e[j]), size: G.e[j] - G.s[j], who: 'Reference', copy: '',
        source: M.gencode ? 'GENCODE ' + M.gencode : 'GENCODE' };
    }
    // What the map draws one by one in view now, by drawsAt(), and the lines of the lanes drawn as a density.
    function inViewRows() {
      var lv = level(), all = lv === 'genome', ch = chrom(), v = view(), rows = [], dense = [];
      function meets(c, s, e) { return all || (c === ch.n && e > v[0] && s < v[1]); }
      if (landed(M.layers.findings || {})) {
        drawn().forEach(function (r) { if (meets(r.chrom, r.s, r.e) && drawsAt('findings', r, lv)) { rows.push(findingRow(r)); } });
      }
      var rl = rohLayer(), banded = 0;
      if (landed(rl)) {
        (all ? M.chroms : [ch]).forEach(function (c) {
          var R = M.roh[c.n];
          for (var i = R && !all ? first(R, v[0]) : 0; R && i < R.n && (all || R.s[i] < v[1]); i++) {
            if (!meets(c.n, R.s[i], R.e[i])) { continue; }
            if (drawsAt('roh', { long: R.o.long[i] }, lv)) { rows.push(rohRow(c, R, i, rl)); } else { banded++; }
          }
        });
      }
      if (banded) {
        dense.push((banded === 1 ? '1 shorter run' : comma(banded) + ' shorter runs') +
          ' of homozygosity, drawn as a faint band: open a chromosome to list them.');
      }
      var gs = !all && drawsAt('genes', null, lv) ? genes[ch.n] : null;
      if (gs && gs.state === 'ready') {
        var G = gs.g, list = [];
        for (var j = first(G, v[0]); j < G.n && G.s[j] < v[1]; j++) { if (G.e[j] > v[0]) { list.push(j); } }
        if (list.length > DENSE) { dense.push(comma(list.length) + ' genes in view, drawn as a density: zoom in to list them.'); }
        else { list.forEach(function (k) { rows.push(geneRow(ch, G, k)); }); }
      }
      return { rows: rows, dense: dense };
    }
    function sortRows(rows) {
      var col = tSort.col, d = tSort.dir;
      function place(a, b) { return a.ci - b.ci || a.s - b.s || a.e - b.e || (a.key < b.key ? -1 : a.key > b.key ? 1 : 0); }
      var by = col === 'where' ? function (a, b) { return (a.layer === 'findings' ? 0 : 1) - (b.layer === 'findings' ? 0 : 1) || place(a, b); }
        : col === 'size' ? function (a, b) { return a.size - b.size || place(a, b); }
          : col === 'kind' ? function (a, b) { return a.rank - b.rank || place(a, b); }
            : function (a, b) { return String(a[col]).localeCompare(String(b[col])) || place(a, b); };
      return rows.sort(function (a, b) { return d * by(a, b); });
    }
    function selected(x) { return x.layer === 'findings' ? sel === x.r.id : x.key === selKey; }
    // On settle: the rows, the count (announced only when it changes) and the density lines.
    function drawTable() {
      var got = inViewRows(), n = got.rows.length;
      tRows = sortRows(got.rows);
      inCount.textContent = n === 1 ? '1 feature' : comma(n) + ' features';
      if (n !== announced) { announced = n; inLive.textContent = inCount.textContent + ' in view'; }
      inDense.textContent = '';
      got.dense.forEach(function (t) { add(inDense, el('p', {}, t)); });
      renderTable();
    }
    function renderTable() {
      // Under 600 px each row is a block of label: value lines, of no fixed height: no scrolled window there,
      // and Show all stops at 200 rows.
      var narrow = !!(window.matchMedia && window.matchMedia('(max-width: 600px)').matches), most = narrow ? T_FULL : T_MAX;
      var n = tRows.length, cap = Math.min(n, most), virtual = showAll && !narrow && cap > T_FULL;
      var sig = [tSort.col, tSort.dir, showAll, narrow, sel, selKey].join('|') + '|' + tRows.map(function (x) { return x.key; }).join('\n');
      if (sig === tSig) { return; }
      tSig = sig;
      var had = document.activeElement, fk = had && inWrap.contains(had) ? had.getAttribute('data-col') ||
        (had.closest('tr') && had.closest('tr').getAttribute('data-key')) : null;
      inWrap.textContent = '';
      inFoot.textContent = '';
      inWrap.classList.toggle('gm-virtual', virtual);
      if (!n) { add(inWrap, el('p', { 'class': 'dim' }, 'No mark drawn one by one in this view.')); return; }
      var tb = el('table', { 'class': 'rows gm-table' }), hr = el('tr', {});
      TCOLS.forEach(function (c) {
        var th = el('th', { scope: 'col' });
        if (tSort.col === c[0]) { th.setAttribute('aria-sort', tSort.dir > 0 ? 'ascending' : 'descending'); }
        var b = el('button', { type: 'button', 'class': 'gm-sort', 'data-col': c[0] }, c[1]);
        b.addEventListener('click', function () {
          if (tSort.col === c[0]) { tSort.dir = -tSort.dir; } else { tSort.col = c[0]; tSort.dir = 1; }
          sortRows(tRows);
          renderTable();
        });
        add(hr, add(th, b));
      });
      var body = el('tbody', {});
      add(tb, add(el('thead', {}), hr), body);
      add(inWrap, tb);
      if (virtual) { windowRows(body, cap); } else {
        for (var i = 0; i < (showAll ? cap : Math.min(n, T_SHORT)); i++) { add(body, rowEl(tRows[i])); }
      }
      if (n > T_SHORT) {
        if (!showAll) {
          add(inFoot, el('span', {}, 'and ' + comma(n - T_SHORT) + ' more in view' + (tSort.col === 'where' && tSort.dir > 0 ?
            ' · sorted by where · findings first' : '')));
        } else if (n > most) {
          add(inFoot, el('span', {}, comma(most) + ' of ' + comma(n) + ' shown; zoom in'));
        }
        add(inFoot, button(showAll ? 'Show fewer' : 'Show all', function () { showAll = !showAll; renderTable(); }, 'gm-showall'));
      }
      var back = fk && (inWrap.querySelector('[data-col="' + fk + '"]') || inWrap.querySelector('tr[data-key="' + fk.replace(/["\\]/g, '\\$&') + '"] button'));
      if (back) { back.focus({ preventScroll: true }); }
    }
    function rowEl(x) {
      var tr = el('tr', { 'data-key': x.key, 'data-layer': x.layer });
      if (selected(x)) { tr.setAttribute('aria-current', 'true'); }
      var cells = [x.layer === 'findings' ? U.findingLevel(x.r.level) : x.kind,
        el('button', { type: 'button', 'class': 'gm-go' }, x.what), x.where, sizeWords(x.size), x.who, x.copy || '—', x.source || '—'];
      cells.forEach(function (c, j) {
        var td = el('td', { 'data-label': TCOLS[j][1] });
        if (j === 2) { td.className = 'mono gm-nw'; } else if (j === 3 || j === 4) { td.className = 'gm-nw'; }
        add(tr, add(td, c));
      });
      return tr;
    }
    // Past 200 rows: only the rows in the scrolled window are in the table, between two spacers.
    function windowRows(body, cap) {
      var top = el('tr', { 'class': 'gm-pad', 'aria-hidden': 'true' }), bot = el('tr', { 'class': 'gm-pad', 'aria-hidden': 'true' });
      add(top, el('td', { colspan: '7' }));
      add(bot, el('td', { colspan: '7' }));
      function fill() {
        var a = clamp(Math.floor(inWrap.scrollTop / ROW_PX) - 10, 0, cap), z = Math.min(cap, a + Math.ceil(inWrap.clientHeight / ROW_PX) + 20);
        body.textContent = '';
        top.firstChild.style.setProperty('height', (a * ROW_PX) + 'px');
        bot.firstChild.style.setProperty('height', ((cap - z) * ROW_PX) + 'px');
        add(body, top);
        for (var i = a; i < z; i++) { add(body, rowEl(tRows[i])); }
        add(body, bot);
      }
      fill();
      inWrap.onscroll = function () { if (!fill.raf) { fill.raf = window.requestAnimationFrame(function () { fill.raf = 0; fill(); }); } };
    }
    // A row flies the map to its feature, a finding with its card: 30% of its length on each side, 20 kb at least.
    function flyRow(key) {
      var x = tRows.filter(function (r) { return r.key === key; })[0];
      if (!x) { return; }
      var len = x.e - x.s, padBp = Math.max(len * 0.3, 20000);
      selKey = x.key;
      openChrom(x.ci, (x.s + x.e) / 2, len + 2 * padBp, true);
      if (x.layer === 'findings') { pick(x.r.id); } else { pick(null); }
    }
    inWrap.addEventListener('click', function (e) {
      var tr = e.target && e.target.closest ? e.target.closest('tbody tr[data-key]') : null;
      if (tr) { flyRow(tr.getAttribute('data-key')); }
    });
    function clearTable() {
      tRows = [];
      tSig = '';
      announced = -1;
      inWrap.textContent = '';
      inDense.textContent = '';
      inFoot.textContent = '';
      inLive.textContent = '';
      inCount.textContent = '';
    }

    // --- loading ---
    function fail(why, words) {
      failure = why;
      // The frame and the margin the last refresh scheduled would draw the old document again: drop them.
      // These lines and the failure guards of draw(), describe() and schedule() cover for each other, so no one
      // of them can go red alone: the failed-reload check undoes them together.
      window.cancelAnimationFrame(raf);
      raf = 0;
      window.clearTimeout(settleT);
      fig.setAttribute('data-pgp-map', why);
      stage.hidden = true;
      marks.hidden = true;
      strip.hidden = true;  // the strip and the table may name a result closed since
      inview.hidden = true;
      clearTable();
      fallback.hidden = false;
      fallback.textContent = '';
      add(fallback, el('p', { 'class': 'gm-fail' }, words));
      var rows = opts.fallback ? opts.fallback() : [];
      if (rows.length) {
        add(fallback, U.table({ columns: ['Level', 'Finding', 'Where'], rows: rows }, 12));
      } else {
        add(fallback, el('p', { 'class': 'dim' }, 'No finding with a place on the map is open on this page; every finding is listed in its section.'));
      }
      narr.textContent = '';
      layersBox.textContent = '';
      chips.textContent = '';  // a chip drawn before the failure may name a result closed since
      marks.textContent = '';  // and so may a mark link, hidden or not
      marksKey = '';
      sel = null;
      card.textContent = '';  // the card, the search's suggestions and Nearest too
      card.hidden = true;
      closeFind();
      nearBtn.hidden = true;
      nearBtn.textContent = '';
      hideTip();
    }
    function load(k) {
      key = k;
      if (typeof window.DecompressionStream === 'undefined' || !window.Response || !window.ReadableStream) {
        fail('nodecompress', 'This browser cannot open the map file. The findings with a place are listed here instead.');
        return;
      }
      fig.setAttribute('data-pgp-map', 'loading');
      fetch(opts.url, { credentials: 'same-origin', cache: 'no-store' }).then(function (r) {
        if (r.status === 401) { window.location.href = '/login?next=' + encodeURIComponent(window.location.pathname); return null; }
        if (!r.ok) {
          return r.text().then(function (t) { throw new Error((t || '').trim() || ('the server answered ' + r.status)); }, function () { throw new Error('the server answered ' + r.status); });
        }
        var t1 = performance.now();
        return new Response(r.body.pipeThrough(new window.DecompressionStream('gzip'))).text().then(function (txt) {
          var p0 = performance.now();
          var d = JSON.parse(txt);
          var m = prepare(d);
          var p1 = performance.now();
          return { m: m, parse: p1 - p0, arrived: t1 };
        });
      }).then(function (x) {
        if (!x || k !== key) { return; }
        var fresh = !M;
        M = x.m;
        failure = null;
        stage.hidden = false;
        marks.hidden = false;
        inview.hidden = false;
        fallback.hidden = true;
        hatch = hatch || pattern(0.16);
        fig.setAttribute('data-pgp-map-parse-ms', String(Math.round(x.parse)));
        fig.removeAttribute('data-pgp-map-paint-ms');
        t0 = x.arrived;
        size();
        if (fresh) {
          tgt.ci = cam.ci = 0;
          tgt.c = cam.c = M.chroms[0].len / 2;
          tgt.ls = cam.ls = Math.log(M.chroms[0].len);
        }
        cacheKey = '';
        drawKey();
        drawPaint();
        drawChips();
        draw();
        describe();
        fig.setAttribute('data-pgp-map', 'drawn');
        if (fresh) { readHash(true); }
      }).catch(function (err) {
        if (k !== key) { return; }
        fail('failed', 'The map could not be loaded: ' + String((err && err.message) || err).split('\n')[0] +
          ' The findings with a place are listed here instead.');
      });
    }

    return {
      node: fig,
      // The page drew again: maybe new results (a new key fetches the document again), maybe a result
      // opened or closed (the veil is asked again).
      refresh: function (k) {
        if (k !== key) { load(k); }
        if (failure) {
          fail(failure, fallback.firstChild ? fallback.firstChild.textContent : '');
          return;
        }
        if (!M) { return; }
        cacheKey = '';
        hideTip();  // a tooltip of a result closed since leaves the page
        drawCard();  // and so does the card of one
        if (!findList.hidden) { runFind(); }  // and a search's suggestion of one
        drawPaint();
        drawChips();
        schedule();
        describe();
      }
    };
  }

  window.PGPMap = { create: create, prepare: prepare, drawsAt: drawsAt, LEVELS: LEVELS };
}());
