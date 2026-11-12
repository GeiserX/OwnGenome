// core.js: what the live site's pages share (home.html, index.html, finding.html, family.html). Element builders,
// a fetch that sends the session and goes to the login on a 401, the pedigree symbol, dates, the
// caveat block, the facts list, the table, and the Journal's level marks. Text only through textContent and createElementNS:
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

  // --- the Journal's shared marks: the level glyph and word ---
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
  // A finding's level: a data-quality finding is always a problem with the data, so its mark is red.
  function findingLevel(lv, words) { return level(lv, words, lv === 'trust' ? 'problem' : undefined); }
  // A plain time: "14:22" in UTC, "4 Nov" for a day.
  var MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  function hm(utc) { var d = new Date(utc); return isFinite(d) ? ('0' + d.getUTCHours()).slice(-2) + ':' + ('0' + d.getUTCMinutes()).slice(-2) : ''; }
  function dm(utc) { var d = new Date(utc); return isFinite(d) ? d.getUTCDate() + ' ' + MONTHS[d.getUTCMonth()] : ''; }
  function hours(a, b) {
    var x = Date.parse(a), y = Date.parse(b);
    return isFinite(x) && isFinite(y) && y > x ? (y - x) / 3600e3 : null;
  }

  // --- the veil: a child's results to act on stay closed until someone opens them on purpose ---------
  // The closed box, the same on a child's page (person.js) and on the family page (family-page.js): what stays
  // closed, and the button that opens it.
  function veil(text, button, onOpen) {
    var v = el('div', { 'class': 'veil' });
    add(v, add(el('div', { 'class': 'vt' }), el('b', {}, 'A result about a child. '), text));
    var b = el('button', { type: 'button', 'class': 'btn' }, button || 'Show the result');
    b.addEventListener('click', onOpen);
    add(v, b);
    return v;
  }
  window.PGPUI = { el: el, svg: svg, add: add, getJSON: getJSON, avatar: avatar, day: day, caveat: caveat,
    facts: facts, table: table, tablist: tablist, firstSentence: firstSentence, level: level, findingLevel: findingLevel,
    LEVEL_WORDS: LEVEL_WORDS, LEVEL_SHORT: LEVEL_SHORT, LEVEL_ORDER: LEVEL_ORDER,
    hm: hm, dm: dm, hours: hours, veil: veil };
}());
