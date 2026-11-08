// core.js: what the live site's pages share (home.html, index.html, finding.html, family.html). Element builders,
// a fetch that sends the session and goes to the login on a 401, the pedigree symbol, dates, the
// caveat block, the facts list, the table, and the Journal's level marks and chromosome map. Text only through textContent and createElementNS:
// no innerHTML, no inline script, no eval. Everything hangs on window.PGPUI.
(function () {
  'use strict';
  var root = document.documentElement;
  var SVGNS = 'http://www.w3.org/2000/svg';
  var blocked = [];
  document.addEventListener('securitypolicyviolation', function (e) {
    blocked.push(e.violatedDirective + ' ' + (e.blockedURI || ''));
    root.setAttribute('data-pgp-blocked', blocked.join(' | '));
  });

  function el(tag, attrs, text) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) {
      if (attrs[k] !== null && attrs[k] !== undefined && attrs[k] !== false) { e.setAttribute(k, attrs[k]); }
    });
    if (text !== undefined && text !== null) { e.textContent = String(text); }
    return e;
  }
  function svg(tag, attrs) {
    var e = document.createElementNS(SVGNS, tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    return e;
  }
  function add(parent) {
    for (var i = 1; i < arguments.length; i++) {
      var c = arguments[i];
      if (c === null || c === undefined || c === false) { continue; }
      parent.appendChild(typeof c === 'string' ? document.createTextNode(c) : c);
    }
    return parent;
  }
  function getJSON(url) {
    return fetch(url, { credentials: 'same-origin', cache: 'no-store' }).then(function (r) {
      if (r.status === 401) { window.location.href = '/login?next=' + encodeURIComponent(window.location.pathname); return null; }
      return r.ok ? r.json() : null;
    });
  }
  // The pedigree symbol: a square (male), a circle (female) or a diamond (not known); filled when sequenced.
  function avatar(sex, filled, ratio) {
    var a = el('span', { 'class': 'avatar sex-' + (sex === 'male' || sex === 'female' ? sex : 'unknown') + (filled ? ' filled' : '') ,
      'aria-hidden': 'true' });
    if (typeof ratio === 'number') {
      a.classList.add('ring');
      a.style.setProperty('--p', String(Math.max(0, Math.min(1, ratio))));
    }
    add(a, el('span', { 'class': 'glyph' }));
    return a;
  }
  function day(utc) { return utc ? String(utc).slice(0, 10) : ''; }
  function firstSentence(t) {
    var m = /^(.+?[.!?])(\s+[A-Z(]|$)/.exec(t || '');
    return m ? m[1] : t;
  }
  // The caveats of a card, word for word: the first sentence of the first one in view, all of them in a
  // details element (in the page, so a reader or a search finds them).
  function caveat(list) {
    if (!list || !list.length) { return null; }
    var box = el('div', { 'class': 'caveat' });
    var lead = firstSentence(list[0]);
    add(box, el('p', { 'class': 'caveat-lead' }, lead));
    var d = el('details', { 'class': 'caveat-full' });
    add(d, el('summary', {}, list.length > 1 || lead !== list[0] ? 'Full caveat' : 'Caveat'));
    list.forEach(function (t) { add(d, el('p', {}, t)); });
    add(box, d);
    return box;
  }
  function facts(rows) {
    if (!rows || !rows.length) { return null; }
    var dl = el('dl', { 'class': 'facts' });
    rows.forEach(function (r) {
      var w = el('div', { 'class': 'fact' });
      add(w, el('dt', {}, r[0]), el('dd', {}, r[1]));
      add(dl, w);
    });
    return dl;
  }
  // A table; on a phone each row becomes a block of label: value pairs (the label is data-label).
  function table(t, limit) {
    if (!t || !t.rows || !t.rows.length) { return null; }
    limit = limit || 8;
    var wrap = el('div', { 'class': 'table-wrap' });
    var tb = el('table', { 'class': 'rows' });
    var head = el('tr', {});
    t.columns.forEach(function (c) { add(head, el('th', { scope: 'col' }, c)); });
    add(tb, add(el('thead', {}), head));
    var body = el('tbody', {});
    t.rows.forEach(function (r, i) {
      var tr = el('tr', i >= limit ? { 'class': 'extra', hidden: 'hidden' } : {});
      r.forEach(function (c, j) { add(tr, el('td', { 'data-label': t.columns[j] || '' }, c)); });
      add(body, tr);
    });
    add(tb, body);
    add(wrap, tb);
    if (t.rows.length > limit) {
      var b = el('button', { type: 'button', 'class': 'more' }, 'Show all ' + t.rows.length + ' rows');
      b.addEventListener('click', function () {
        Array.prototype.forEach.call(body.querySelectorAll('tr.extra'), function (x) { x.hidden = false; });
        b.remove();
      });
      add(wrap, b);
    }
    // Rows past results.json's cap are only in the report: said apart from the button, never promised by it.
    if (t.more > 0) {
      add(wrap, el('p', { 'class': 'dim' }, t.more + ' more rows in the report'));
    }
    return wrap;
  }
  // role=tablist keyboard pattern: arrow keys move, Home and End jump; Enter or a click activates.
  function tablist(list, onPick) {
    list.addEventListener('keydown', function (e) {
      var tabs = Array.prototype.slice.call(list.querySelectorAll('[role=tab]'));
      var i = tabs.indexOf(document.activeElement);
      if (i < 0) { return; }
      var j = { ArrowDown: i + 1, ArrowRight: i + 1, ArrowUp: i - 1, ArrowLeft: i - 1, Home: 0, End: tabs.length - 1 }[e.key];
      if (j === undefined) { return; }
      e.preventDefault();
      j = (j + tabs.length) % tabs.length;
      tabs[j].focus();
      onPick(tabs[j]);
    });
  }

  // --- the Journal's shared marks: the level glyph and word, and the chromosome map of findings ---
  var LEVEL_WORDS = { act: 'Take to a doctor', medicine: 'Tell your prescriber', review: 'Worth a look',
    carrier: 'For family planning', trust: 'Data quality', info: 'For information' };
  var LEVEL_SHORT = { act: 'doctor', medicine: 'prescriber', review: 'look', carrier: 'family', trust: 'quality', info: 'info' };
  var LEVEL_ORDER = ['act', 'medicine', 'review', 'carrier', 'trust', 'info'];
  // A level is always its glyph and a word, never a colour alone.
  function level(lv, words, extra) {
    var s = el('span', { 'class': 'lv ' + lv + (extra ? ' ' + extra : ''), 'data-level': lv });
    add(s, el('span', { 'class': 'g', 'aria-hidden': 'true' }), words === undefined ? LEVEL_WORDS[lv] || lv : words);
    return s;
  }
  // GRCh38: [name, length Mb, centromere Mb].
  var CHR = [['1', 248.96, 123.4], ['2', 242.19, 93.9], ['3', 198.3, 90.9], ['4', 190.21, 50], ['5', 181.54, 48.8],
    ['6', 170.81, 59.8], ['7', 159.35, 60.1], ['8', 145.14, 45.2], ['9', 138.39, 43], ['10', 133.8, 39.8],
    ['11', 135.09, 53.4], ['12', 133.28, 35.5], ['13', 114.36, 17.7], ['14', 107.04, 17.2], ['15', 101.99, 19],
    ['16', 90.34, 36.8], ['17', 83.26, 25.1], ['18', 80.37, 18.5], ['19', 58.62, 26.2], ['20', 64.44, 28.1],
    ['21', 46.71, 12], ['22', 50.82, 15], ['X', 156.04, 60.6], ['Y', 57.23, 10.4]];
  var LV_COLOUR = { act: '#d03b3b', medicine: '#ec835a', review: '#fab219', carrier: '#9085e9', trust: '#c3c2b7', info: '#c3c2b7' };
  function glyph(g, lv, x, y) {
    var c = LV_COLOUR[lv] || LV_COLOUR.info;
    if (lv === 'act') { g.appendChild(svg('rect', { x: x - 4, y: y - 4, width: 8, height: 8, fill: c, transform: 'rotate(45 ' + x + ' ' + y + ')' })); }
    else if (lv === 'medicine') { g.appendChild(svg('rect', { x: x - 5.5, y: y - 3, width: 11, height: 6, rx: 3, fill: c })); }
    else if (lv === 'review') { g.appendChild(svg('path', { d: 'M' + x + ' ' + (y - 5) + 'L' + (x + 5) + ' ' + (y + 4) + 'L' + (x - 5) + ' ' + (y + 4) + 'Z', fill: c })); }
    else if (lv === 'carrier') {
      g.appendChild(svg('circle', { cx: x - 2.6, cy: y, r: 3.4, fill: 'none', stroke: c, 'stroke-width': 1.6 }));
      g.appendChild(svg('circle', { cx: x + 2.6, cy: y, r: 3.4, fill: 'none', stroke: c, 'stroke-width': 1.6 }));
    } else { g.appendChild(svg('circle', { cx: x, cy: y, r: 4, fill: 'none', stroke: c, 'stroke-width': 1.5 })); }
  }
  // "chr16:2087914" or "chr16:2084000-2096000" -> {chr: '16', mb: 2.08}; null for anything else.
  function locus(v) {
    var m = /^(?:chr)?([0-9]{1,2}|X|Y):([0-9]+)/.exec(String(v || ''));
    return m ? { chr: m[1], mb: Number(m[2]) / 1e6 } : null;
  }
  // The chromosome map: the 24 chromosomes in a row, each finding a glyph on its chromosome, with its
  // name beside it and the same words in a title (and in the list under the figure).
  function karyogram(marks, opts) {
    opts = opts || {};
    var W = opts.width || 680, H = opts.height || 190, top = 22, bottom = 22;
    var slot = (W - 56) / CHR.length, scale = (H - top - bottom) / 249;
    var s = svg('svg', { viewBox: '0 0 ' + W + ' ' + H, role: 'img', 'aria-label': opts.aria || 'Where the findings sit, by chromosome' });
    CHR.forEach(function (c, i) {
      var x = i * slot + slot / 2 - 3.5;
      s.appendChild(svg('rect', { 'class': 'chr', x: x, y: top, width: 7, height: Math.max(1, c[2] * scale - 1.2), rx: 3.5 }));
      s.appendChild(svg('rect', { 'class': 'chr', x: x, y: top + c[2] * scale + 1.2, width: 7, height: Math.max(1, (c[1] - c[2]) * scale - 1.2), rx: 3.5 }));
      var t = svg('text', { 'class': 'lab', x: x + 3.5, y: H - 6, 'text-anchor': 'middle' });
      t.textContent = c[0];
      s.appendChild(t);
    });
    (marks || []).forEach(function (m) {
      var i = -1;
      CHR.forEach(function (c, j) { if (c[0] === m.chr) { i = j; } });
      if (i < 0) { return; }
      var x = i * slot + slot / 2 + 7, y = top + Math.min(m.mb, CHR[i][1]) * scale;
      var g = svg('g', { 'class': 'mk', tabindex: '0' });
      var tt = svg('title', {});
      tt.textContent = (LEVEL_WORDS[m.lv] || '') + ': ' + (m.tip || m.label);
      g.appendChild(tt);
      g.appendChild(svg('line', { x1: x - 10.5, x2: x - 1, y1: y, y2: y, stroke: LV_COLOUR[m.lv] || LV_COLOUR.info, 'stroke-width': 2 }));
      glyph(g, m.lv, x + 4, y);
      if (m.label) {
        var t = svg('text', { x: x + 11, y: y + 3.5, stroke: '#0e0f10', 'stroke-width': 4, 'paint-order': 'stroke', 'stroke-linejoin': 'round' });
        t.textContent = m.label;
        g.appendChild(t);
      }
      if (m.href) {
        g.setAttribute('role', 'link');
        g.addEventListener('click', function () { window.location.href = m.href; });
        g.addEventListener('keydown', function (e) { if (e.key === 'Enter') { window.location.href = m.href; } });
      }
      s.appendChild(g);
    });
    return s;
  }
  // A plain time: "14:22" in UTC, "4 Nov" for a day.
  var MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  function hm(utc) { var d = new Date(utc); return isFinite(d) ? ('0' + d.getUTCHours()).slice(-2) + ':' + ('0' + d.getUTCMinutes()).slice(-2) : ''; }
  function dm(utc) { var d = new Date(utc); return isFinite(d) ? d.getUTCDate() + ' ' + MONTHS[d.getUTCMonth()] : ''; }
  function hours(a, b) {
    var x = Date.parse(a), y = Date.parse(b);
    return isFinite(x) && isFinite(y) && y > x ? (y - x) / 3600e3 : null;
  }

  window.PGPUI = { el: el, svg: svg, add: add, getJSON: getJSON, avatar: avatar, day: day, caveat: caveat,
    facts: facts, table: table, tablist: tablist, firstSentence: firstSentence, level: level, glyph: glyph,
    karyogram: karyogram, locus: locus, LEVEL_WORDS: LEVEL_WORDS, LEVEL_SHORT: LEVEL_SHORT, LEVEL_ORDER: LEVEL_ORDER,
    hm: hm, dm: dm, hours: hours };
}());
