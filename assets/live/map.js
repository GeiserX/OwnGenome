// map.js: Figure 1 of a person's page, the genome map. It draws the map document the server builds with
// `og map build` (GET map.json beside the page: gzip-wrapped columnar JSON, rust/crates/og-map/README.md):
// the 24 chromosomes as rows with their bands, the masked spans hatched, each finding with a place as its
// level's glyph, and the runs of homozygosity; a click opens a chromosome in place, where the bands, the
// findings and the runs get lanes of their own under a ruler. Every mark comes from the document, never from
// results.json, and every sentence the margin shows is one the document wrote (text, phrase, empty, note),
// with the counts filled in here.
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
  // The chromosome level ends at 16,667 bp per CSS pixel: closer views (region, gene, reads) come later.
  var CHROM_BPP = 16667;
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
        gene: f.gene[i], title: f.title[i], href: f.href[i], phrase: f.phrase[i], veil: f.veil[i] >= 0 ? veilIds[f.veil[i]] : null });
    });
    var gated = {};
    rows.forEach(function (r) { if (r.gate) { gated[r.id] = r; } });
    return { doc: d, id: id, label: person.label || id, chroms: chroms, by: by, bands: bands, roh: roh, bins: bins,
      binSize: (per.roh_bins && per.roh_bins.size) || 1e6, findings: rows, gated: gated, veilIds: veilIds,
      layers: d.layers || {}, text: d.text || {}, art: (d.dict && d.dict.art) || [] };
  }

  // --- the figure ---------------------------------------------------------------------------------------
  // create({url, veil(id, gated) -> 'open'|'closed'|'opened', now(id) -> shown at this moment, open(id),
  // close(id), fallback() -> rows}): one figure per page. The page calls refresh(key) on each of its redraws; a new key
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
    add(gauge, stopG, stopC);
    var crumb = el('nav', { 'class': 'gm-crumb mono', 'aria-label': 'Where the map is' });
    var chips = el('div', { 'class': 'gm-chips' });
    add(head, kicker, add(el('div', { 'class': 'gm-titlerow' }), title, gauge), crumb, chips);
    var stage = el('div', { 'class': 'gm-stage' });
    var canvas = el('canvas', { 'class': 'gm-canvas', tabindex: '0', role: 'img', 'aria-describedby': 'gm-narration' });
    var tip = el('div', { 'class': 'gm-tip', role: 'tooltip', hidden: 'hidden' });
    var toast = el('div', { 'class': 'gm-toast', role: 'status', hidden: 'hidden' });
    add(stage, canvas, tip, toast);
    // The marks drawn, as links: the keyboard and screen-reader path to each finding (Tab steps through them).
    var marks = el('ul', { 'class': 'gm-marks', 'aria-label': 'Findings on the map' });
    var fallback = el('div', { 'class': 'gm-fallback', hidden: 'hidden' });
    add(main, head, stage, marks, fallback);
    var margin = el('aside', { 'class': 'gm-margin' });
    var mkick = el('div', { 'class': 'fig-label' }, 'The whole genome');
    var mline = el('div', { 'class': 'gm-mline mono' });
    var narr = el('div', { 'class': 'gm-narration', id: 'gm-narration', 'aria-live': 'polite' });
    var layersBox = el('dl', { 'class': 'gm-layers' });
    var keyBox = el('div', { 'class': 'gm-key', 'aria-label': 'Key' });
    // The layers and the key: open beside the map, folded under it on a phone.
    var more = el('details', { 'class': 'gm-more' });
    more.open = !(window.matchMedia && window.matchMedia('(max-width: 700px)').matches);
    add(more, el('summary', {}, 'Layers and key'), layersBox, keyBox);
    add(margin, mkick, mline, narr, more);
    add(fig, main, margin);

    var ctx = canvas.getContext('2d');
    var CW = 0, CH = 0, DPR = 1, phone = false;
    var cam = { ci: 0, c: 0, ls: 0, k: 0 }, tgt = { ci: 0, c: 0, ls: 0, k: 0 };
    var hover = { row: -1, hit: null }, hits = [], raf = 0, settleT = 0, hashT = 0, confirming = null, touched = false;
    var cache = null, cacheKey = '';
    var hatch = null;

    // --- the veil: what this page keeps closed at this moment ---
    function stateOf(id) { return opts.veil(id, !!(M && M.gated[id])); }
    // A row shows when its result is not closed and the page shows its finding at this moment (a replay
    // before the finding's step ended does not).
    function shownRow(r) { return (!r.veil || stateOf(r.veil) !== 'closed') && (!opts.now || opts.now(r.id)); }
    function drawn() { return M ? M.findings.filter(shownRow) : []; }

    // --- geometry ---
    function chrom() { return M.chroms[cam.ci]; }
    function rowH() { var lo = phone ? ROW_PHONE : ROW_MIN; return M ? clamp((CH - 18) / M.chroms.length, lo, 44) : lo; }
    function gutG() { return phone ? 30 : 40; }
    function pad() { return phone ? 10 : 22; }
    function gut() { return lerp(gutG(), phone ? 64 : 128, cam.k); }
    function plotW() { return CW - gut() - pad(); }
    function minSpan(len) { return Math.min(len, CHROM_BPP * Math.max(200, CW - (phone ? 64 : 128) - pad())); }
    function view() {
      var span = Math.exp(cam.ls);
      return [cam.c - span / 2, cam.c + span / 2, span];
    }
    function gx(bp) { return gutG() + bp * (CW - gutG() - pad()) / M.chroms[0].len; }
    function ox(bp) { var v = view(); return gut() + (bp - v[0]) / v[2] * plotW(); }
    function xAt(bp, open) { return open ? lerp(gx(bp), ox(bp), cam.k) : gx(bp); }
    function rowY(i) { return 9 + i * rowH() + rowH() / 2; }
    function level() {
      if (cam.k < 0.5) { return 'genome'; }
      return 'chromosome';
    }
    // The lanes of an open chromosome, top to bottom: y of each, stretched to the canvas.
    function lanes() {
      var free = Math.max(0, CH - 230), st = 1 + Math.min(0.75, free / 400);
      var fy = 118, fh = 50 * Math.min(1.2, st), ry = fy + fh + 14, rh = 22 * st;
      return { over: 10, ruler: 46, bands: 72, bh: 38, fy: fy, fh: fh, ry: ry, rh: rh };
    }

    // --- the canvas size, the hatch pattern ---
    // False while the stage has no width (not in the page yet, or hidden): sized then, a phone would get the
    // desktop's rows. The ResizeObserver draws once it has one.
    function size() {
      if (!stage.clientWidth) { return false; }
      var w = Math.max(280, Math.floor(stage.clientWidth));
      phone = w < 520;
      var n = M ? M.chroms.length : 24;
      var h = 18 + n * (phone ? ROW_PHONE : ROW_MIN);
      DPR = Math.min(2, window.devicePixelRatio || 1);
      if (w !== CW || h !== CH) {
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
      var fs = drawn();
      var occupied = {};
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
        var bins = M.bins[ch.n], ry = yc + 6.5;
        if (bins) {
          c.fillStyle = INK2;
          for (var b = 0; b < bins.length; b++) {
            if (!bins[b]) { continue; }
            c.globalAlpha = 0.22 * Math.min(1, bins[b] / M.binSize);
            c.fillRect(gx(b * M.binSize), ry - 1, Math.max(1, gx((b + 1) * M.binSize) - gx(b * M.binSize)), 2);
          }
          c.globalAlpha = 1;
        }
        var rr = M.roh[ch.n];
        if (rr) {
          c.strokeStyle = INK;
          c.lineWidth = 2.5;
          for (var j = 0; j < rr.n; j++) {
            if (!rr.o.long[j]) { continue; }
            c.beginPath(); c.moveTo(gx(rr.s[j]), ry); c.lineTo(gx(rr.e[j]), ry); c.stroke();
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
      if (!M) { return; }
      var t = performance.now();
      if (!size()) { return; }  // the page may have changed width since the last frame (a phone turned, a fold opened)
      step();
      ctx.setTransform(1, 0, 0, 1, 0, 0);
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
      hits = [];
      if (cam.k < 0.999) { drawGenome(); }
      if (cam.k > 0.001) { drawOpen(); }
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
        M.chroms.forEach(function (c, i) {
          hits.push({ kind: 'row', i: i, x0: 0, x1: CW, y0: rowY(i) - rowH() / 2, y1: rowY(i) + rowH() / 2 });
        });
        drawn().forEach(function (r) {
          hits.push({ kind: 'finding', r: r, x: gx((r.s + r.e) / 2), y: rowY(M.by[r.chrom]) - 1, rad: 7 });
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
        ctx.fillText(mb(at, stp < 1e6 ? 1 : 0) + ' Mb', x, L.ruler + 8);
      }
      // Band names where they fit.
      var b = M.bands[ch.n];
      if (b) {
        ctx.font = '10px ' + SANS;
        ctx.textBaseline = 'middle';
        for (var i = first(b, v[0]); i < b.n && b.s[i] < v[1]; i++) {
          var a = ox(b.s[i]), z = ox(b.e[i]), nm = b.o.name[i];
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
      ctx.textAlign = 'left';
      ctx.save();
      ctx.beginPath(); ctx.rect(g - 6, 0, W + 12, CH); ctx.clip();
      openFindings(L, v, ch);
      openRoh(L, v, ch);
      // The reading line: a dashed brass line through every lane at the centre, its place in a pill.
      var cx = g + W / 2;
      ctx.setLineDash([3, 3]);
      ctx.strokeStyle = BRASS;
      ctx.globalAlpha = k * 0.8;
      ctx.beginPath(); ctx.moveTo(cx, L.ruler + 22); ctx.lineTo(cx, L.ry + L.rh + 6); ctx.stroke();
      ctx.setLineDash([]);
      ctx.restore();
      var pill = mb(cam.c, 1) + ' Mb', pw = measure('600 10.5px ' + SANS, pill) + 12;
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
      var inView = fs.filter(function (r) { return r.e > v[0] && r.s < v[1]; });
      if (!inView.length) {
        laneWords(ctx, fs.length ? 'No finding in this stretch.' : noFindingHere(), g + 4, L.fy + 16);
        return;
      }
      var used = [];
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
        ctx.globalAlpha = cam.k;
        glyph(ctx, r.level, x, y, 6);
        hits.push({ kind: 'finding', r: r, x: x, y: y, rad: 9 });
        var name = r.gene || r.title, words = U.LEVEL_WORDS[r.level] || r.level;
        var w = Math.max(measure('500 12px ' + SANS, name), measure('10.5px ' + SANS, words));
        var lx = x + 11;
        if (lx + w > CW - 4 || used.some(function (u) { return !(lx + w < u[0] || lx > u[1]); })) { return; }
        used.push([lx - 4, lx + w + 6]);
        ctx.textBaseline = 'middle';
        label(ctx, name, lx, y - 1, '500 12px ' + SANS, INK);
        label(ctx, words, lx, y + 13, '10.5px ' + SANS, INK3);
      });
    }
    function openRoh(L, v, ch) {
      var lay = M.layers.roh || {}, g = gut(), y = L.ry + L.rh / 2;
      if (!landed(lay)) { laneWords(ctx, layerWords(lay), g + 4, y); return; }
      var r = M.roh[ch.n];
      if (!r || !r.n) { laneWords(ctx, 'No run of homozygosity on this chromosome.', g + 4, y); return; }
      for (var i = first(r, v[0]); i < r.n && r.s[i] < v[1]; i++) {
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
    function noFindingHere() {
      return closedAny() ? 'No open finding with a place on this chromosome.' : 'No finding with a place on this chromosome.';
    }
    function closedAny() { return M.veilIds.some(function (id) { return stateOf(id) === 'closed'; }); }
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
    function schedule() { if (!raf) { raf = window.requestAnimationFrame(function () { last = last || performance.now(); draw(); }); } }
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
      return ch.n + ':' + mb(tgt.c - span / 2) + '-' + mb(tgt.c + span / 2) + 'M';
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
      if (!h) { return; }
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
      settleT = window.setTimeout(describe, 140);
    }
    function describe() {
      if (!M) { return; }
      var lv = level(), fs = drawn(), ch = chrom();
      fig.setAttribute('data-level', lv);
      fig.setAttribute('data-chrom', lv === 'genome' ? '' : ch.n);
      fig.setAttribute('data-view', locusWords());
      stopG.setAttribute('aria-pressed', String(lv === 'genome'));
      stopC.setAttribute('aria-pressed', String(lv !== 'genome'));
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
        title.textContent = 'Chromosome ' + short(ch.n);
        kicker.textContent = 'Figure 1 · Genome map · ' + M.label;
        add(crumb, el('span', { 'class': 'sep', 'aria-hidden': 'true' }, ' › '), el('span', { 'aria-current': 'location' }, ch.n));
        mkick.textContent = 'Under the line';
        mline.textContent = ch.n + ' · ' + mb(cam.c, 2) + ' Mb' + bandAt(ch, cam.c);
      }
      canvas.setAttribute('aria-label', lv === 'genome' ? 'The genome map: ' + nY + ' chromosomes with ' + fs.length + ' findings drawn'
        : 'Chromosome ' + short(ch.n) + ' open on the genome map');
      narr.textContent = '';
      var v = view();
      var shown = lv === 'genome' ? fs : fs.filter(function (r) { return r.chrom === ch.n && r.e > v[0] && r.s < v[1]; });
      var forms = (M.text.findings || {});
      if (lv === 'genome') {
        var form = fs.length === 0 ? forms.none : fs.length === 1 ? forms.one : forms.many;
        add(narr, runs((form || []).map(function (x) { return { t: x.t.replace('{n}', String(fs.length)), b: x.b }; })));
      } else if (M.text.chrom && M.text.chrom[ch.n]) {
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
        var tail = runs(M.text.genome), lead = narr.children.length === 1 ? narr.firstChild : null;
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
        if (n && landed(M.layers.roh || {})) {
          add(narr, el('p', {}, 'Runs of homozygosity in view: ' + n + ', the longest ' + mb(longest) + ' Mb.'));
        }
        if (!phone) { add(narr, el('p', { 'class': 'dim' }, 'Drag to move along it, scroll to zoom, Esc for the whole genome.')); }
      }
      drawLayers();
      drawMarks();
    }
    function bandAt(ch, bp) {
      var b = M.bands[ch.n];
      if (!b) { return ''; }
      for (var i = first(b, bp); i < b.n && b.s[i] <= bp; i++) {
        if (b.e[i] > bp) { return ' · band ' + short(ch.n) + b.o.name[i]; }
      }
      return '';
    }
    function drawLayers() {
      layersBox.textContent = '';
      [['findings', 'Findings'], ['roh', 'Runs of homozygosity']].forEach(function (x) {
        var lay = M.layers[x[0]];
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
    function drawMarks() {
      marks.textContent = '';
      var lv = level(), ch = chrom(), v = view();
      drawn().filter(function (r) { return lv === 'genome' || (r.chrom === ch.n && r.e > v[0] && r.s < v[1]); })
        .sort(function (a, b) { return M.by[a.chrom] - M.by[b.chrom] || a.s - b.s; }).forEach(function (r) {
          var a = el('a', { href: r.href, 'data-id': r.id, 'data-level': r.level, 'data-chrom': r.chrom,
            'data-s': String(r.s), 'data-e': String(r.e) }, (U.LEVEL_WORDS[r.level] || r.level) + ': ' + r.title + ', ' + r.chrom);
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

    // --- the closed results: one dashed chip per kind, and the opened ones ---
    function drawChips() {
      chips.textContent = '';
      if (!M) { return; }
      var acts = [], gates = [], opened = [];
      M.veilIds.forEach(function (id) {
        var s = stateOf(id);
        if (s === 'closed') { (M.gated[id] ? gates : acts).push(id); } else if (s === 'opened') { opened.push(id); }
      });
      // Why it is closed: said beside the count, left out on a phone (the card says it too).
      function why() { return el('span', { 'class': 'gm-why' }, ' · Not drawn until opened on purpose'); }
      if (acts.length) {
        var c = el('div', { 'class': 'gm-chip closed', 'data-closed': 'act' });
        add(c, add(el('span', {}, acts.length === 1 ? '1 result closed' : acts.length + ' results closed'), why()));
        add(c, confirming === 'act' ? confirmBox(acts[0], false) : button('Open', function () { confirming = 'act'; drawChips(); }));
        add(chips, c);
      }
      gates.forEach(function (id) {
        var g = M.gated[id], name = g.gene || 'This result';
        var c = el('div', { 'class': 'gm-chip closed', 'data-closed': 'gate' });
        add(c, add(el('span', {}, name + ' closed'), why()));
        add(c, confirming === id ? confirmBox(id, true) : button('Open ' + name, function () { confirming = id; drawChips(); }));
        add(chips, c);
      });
      opened.forEach(function (id) {
        var r = M.findings.filter(function (x) { return x.id === id; })[0];
        var c = el('div', { 'class': 'gm-chip opened', 'data-opened': '1' });
        add(c, el('span', {}, 'Opened: ' + (r ? r.title : 'a result to act on')));
        add(c, button('Close it again', function () { confirming = null; opts.close(id); }));
        add(chips, c);
      });
    }
    function confirmBox(id, gated) {
      var box = el('span', { 'class': 'gm-confirm' });
      add(box, el('span', {}, gated ? 'It is kept behind its warning on every page. Open it here and in its section?'
        : 'A result to act on in adult life is usually not looked at before the child can choose. Open it here and in its section?'));
      add(box, button('Show the result', function () { confirming = null; opts.open(id); }, 'gm-yes'));
      add(box, button('Keep it closed', function () { confirming = null; drawChips(); }, 'gm-no'));
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
      if (h.kind === 'finding') { return [h.r.phrase, 'Click for the finding’s page']; }
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
      canvas.style.setProperty('cursor', h && (h.kind === 'finding' || h.kind === 'row') ? 'pointer' : level() === 'genome' ? 'default' : 'grab');
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
      if (h && h.kind === 'finding') { window.location.href = h.r.href; return; }
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
    canvas.addEventListener('keydown', function (e) {
      if (!M || e.altKey || e.ctrlKey || e.metaKey) { return; }
      var open = level() !== 'genome', span = Math.exp(tgt.ls);
      var done = true;
      if (e.key === 'Escape' || e.key === '0') { if (open) { closeChrom(true); } else { done = false; } }
      else if (e.key === '[' || e.key === ']') {
        var i = (open ? tgt.ci : (e.key === ']' ? -1 : M.chroms.length)) + (e.key === ']' ? 1 : -1);
        openChrom((i + M.chroms.length) % M.chroms.length, undefined, undefined, true);
      } else if (open && (e.key === '+' || e.key === '=' || e.key === '-' || e.key === '−')) {
        tgt.ls = Math.log(span * (e.key === '-' || e.key === '−' ? 1.6 : 1 / 1.6));
        clampTgt();
        changed(false);
      } else if (open && (e.key === 'ArrowLeft' || e.key === 'ArrowRight')) {
        tgt.c += (e.key === 'ArrowLeft' ? -0.2 : 0.2) * span;
        clampTgt();
        changed(false);
      } else if (!open && e.key === 'Enter') { openChrom(tgt.ci, undefined, undefined, true); }
      else { done = false; }
      if (done) { e.preventDefault(); }
    });
    stopG.addEventListener('click', function () { closeChrom(true); });
    stopC.addEventListener('click', function () {
      if (level() !== 'genome') { openChrom(tgt.ci, undefined, undefined, true); return; }
      // From the genome: the chromosome of the first finding by level, else chromosome 1.
      var f = drawn().slice().sort(function (a, b) { return PRI[a.level] - PRI[b.level] || M.by[a.chrom] - M.by[b.chrom]; })[0];
      openChrom(f ? M.by[f.chrom] : 0, undefined, undefined, true);
    });
    window.addEventListener('hashchange', function () { readHash(false); });
    window.addEventListener('popstate', function () { readHash(false); });
    if (window.ResizeObserver) { new window.ResizeObserver(function () { if (M) { size(); clampTgt(); schedule(); } }).observe(stage); }

    // --- loading ---
    function fail(why, words) {
      failure = why;
      fig.setAttribute('data-pgp-map', why);
      stage.hidden = true;
      marks.hidden = true;
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
        drawChips();
        schedule();
        describe();
      }
    };
  }

  window.PGPMap = { create: create, prepare: prepare };
}());
