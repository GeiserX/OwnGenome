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

  window.PGPUI = { el: el, svg: svg, add: add, getJSON: getJSON, avatar: avatar, day: day, caveat: caveat,
    facts: facts, table: table, tablist: tablist, firstSentence: firstSentence, level: level, findingLevel: findingLevel,
    LEVEL_WORDS: LEVEL_WORDS, LEVEL_SHORT: LEVEL_SHORT, LEVEL_ORDER: LEVEL_ORDER,
    hm: hm, dm: dm, hours: hours };
  // --- the report gate ----------------------------------------------------------------------------------
  // A child's report stays closed while any of the child's results to act on (findings of level "act") is
  // unopened. The page is then a function of the opened results alone and of how many are closed: no view (person.js,
  // finding.js, map.js, family-page.js, home.js, live.js, the plugins) receives the child's report, because every
  // view reads through fetch and EventSource and this gate holds both. Adults and children with no result to act on
  // pass through untouched: each file is fetched as often, and in the same order, as without the gate.
  var originalFetch = window.fetch.bind(window);
  var records = {}, current = null, ownsPage = false, eventSources = [], shared = {}, cleaned = false;
  var CLOSED_TEXT = "This child's report stays closed until every result to act on in adult life is opened on purpose.";
  var TRANSFER = 'og-report-opening', NONCE = 'og-report-opening:';
  var findingPath = /\/finding\//.test(location.pathname) ? location.pathname : null;
  // The opening ticket of the navigation that replaced the last document (transfer below): read and removed at
  // once, before any network load, and valid only for this navigation (its one-use nonce in window.name).
  var ticket = null;
  try {
    var text = sessionStorage.getItem(TRANSFER), t = text ? JSON.parse(text) : null;
    sessionStorage.removeItem(TRANSFER);
    if (window.name.indexOf(NONCE) === 0) {
      if (t && typeof t.nonce === 'string' && window.name === NONCE + t.nonce && t.path === location.pathname) { ticket = t; }
      window.name = t && typeof t.name === 'string' ? t.name : '';
    }
  } catch (e) { ticket = null; }
  function absolute(url, base) { return new URL(url, base || location.href).href; }
  function copy(x) { return x === undefined ? null : JSON.parse(JSON.stringify(x)); }
  function read(r) {
    if (r.status === 401) { window.location.href = '/login?next=' + encodeURIComponent(location.pathname); }
    return r.ok ? r.json() : null;
  }
  function jsonAt(url) {
    return originalFetch(url, { credentials: 'same-origin', cache: 'no-store' }).then(read).catch(function () { return null; });
  }
  // A file the gate reads before the views: kept once, so the view's own first read of it costs no request.
  function sharedJSON(url) {
    url = absolute(url);
    return originalFetch(url, { credentials: 'same-origin', cache: 'no-store' }).then(function (r) {
      return r.text().then(function (body) {
        if (r.ok) { shared[url] = { body: body, status: r.status, type: r.headers.get('Content-Type') || 'application/json' }; }
        return read({ status: r.status, ok: r.ok, json: function () { return Promise.resolve(JSON.parse(body)); } });
      });
    }).catch(function () { return null; });
  }
  function random() { var b = new Uint32Array(1); crypto.getRandomValues(b); return b[0] / 4294967296; }
  // The results to act on, one record each: the findings of level "act" (bin/live_results.py). null when unreadable.
  function actsOf(r) {
    if (!r || !Array.isArray(r.findings)) { return null; }
    var acts = [], seen = {};
    r.findings.forEach(function (f) {
      if (f && f.level === 'act' && typeof f.id === 'string' && !seen[f.id]) { seen[f.id] = true; acts.push(f); }
    });
    return acts;
  }
  function byId(rec, id) { return (rec.acts || []).filter(function (f) { return f.id === id; })[0]; }
  function same(a, b) { return !!a && !!b && JSON.stringify(a) === JSON.stringify(b); }
  // perm: the acts in a random order drawn for this page load. A closed control is named by its place in it
  // (closed-1, closed-2), so no control's place says anything about the act behind it. order: the opened acts in
  // the order they were opened (opened-1, opened-2).
  function shuffled(ids) {
    var out = ids.slice();
    for (var i = out.length - 1; i > 0; i--) { var j = Math.floor(random() * (i + 1)), x = out[i]; out[i] = out[j]; out[j] = x; }
    return out;
  }
  function makeRecord(url, child, r, old) {
    var acts = actsOf(r), rec = { url: url, child: child, acts: acts, opened: {}, order: [], perm: [] };
    if (!acts) { return rec; }
    // An act whose record changed is closed again; an unchanged one keeps its opening for this page.
    rec.order = old ? old.order.filter(function (id) { return same(byId(old, id), byId(rec, id)); }) : [];
    rec.order.forEach(function (id) { rec.opened[id] = true; });
    var ids = acts.map(function (f) { return f.id; });
    rec.perm = old && same(ids, (old.acts || []).map(function (f) { return f.id; })) ? old.perm : shuffled(ids);
    return rec;
  }
  function record(url, child, cached) {
    url = absolute(url);
    return (cached ? sharedJSON : jsonAt)(url + 'results.json').then(function (r) {
      var old = records[url];
      var rec = makeRecord(url, child !== undefined ? child : old ? old.child : true, r, old);
      records[url] = rec;
      return rec;
    });
  }
  function shut(rec) { return rec.child && (!rec.acts || rec.acts.some(function (f) { return !rec.opened[f.id]; })); }
  function projection(rec) {
    var closed = rec.perm.filter(function (id) { return !rec.opened[id]; });
    return { closed: closed.map(function (id, k) { return { handle: 'closed-' + (k + 1), id: id }; }),
      opened: rec.order.map(function (id, k) { return { handle: 'opened-' + (k + 1), id: id, finding: copy(byId(rec, id)) }; }) };
  }
  function handles(rec) { var p = projection(rec); return p.closed.concat(p.opened); }
  // The address of a closed page names a handle or nothing: any other #f/<x> (an id, a guess) and any query
  // land on the first closed handle, so a typed address tells nothing about which result is closed.
  function normalizeHash(rec) {
    if (!location.hash && !location.search) { return; }
    var all = handles(rec), h = location.hash.slice(3);
    var found = /^#f\/(closed|opened)-[1-9][0-9]*$/.test(location.hash) && all.some(function (a) { return a.handle === h; });
    var to = found ? location.hash : all.length ? '#f/' + all[0].handle : '';
    if (location.hash !== to || location.search) { history.replaceState(null, '', location.pathname + to); }
  }
  function actRecord(f) {
    var box = el('div', { 'class': 'opened-act' });
    add(box, el('h2', {}, f.title), el('p', {}, f.detail), el('p', {}, f.text));
    var plan = el('ol', {});
    (f.plan || []).forEach(function (t) { add(plan, el('li', {}, t)); });
    add(box, plan, facts(f.evidence), caveat(f.caveat));
    return box;
  }
  // The controls of a closed report: the fixed text, the count, one block of identical closed controls, then each
  // opened act's own record, in the order it was opened. PREFIX keeps element ids apart on a page with two children.
  function controls(rec, changed, prefix) {
    var p = projection(rec), box = el('section', { 'class': 'report-gate' }), n = (rec.acts || []).length;
    add(box, el('p', { 'class': 'lede' }, CLOSED_TEXT));
    add(box, el('p', {}, rec.acts ? n + (n === 1 ? ' result' : ' results') + ' to act on' : 'The results cannot be read yet. Reload to try again.'));
    function row(a, open) {
      var art = el('article', { id: prefix + 'f-' + a.handle, 'data-finding': a.handle, 'class': 'finding' });
      if (open) { add(art, actRecord(a.finding)); }
      var b = el('button', { type: 'button', 'class': 'btn', 'data-opening': a.handle }, open ? 'Close it again' : 'Show the result');
      b.addEventListener('click', function () { toggle(rec, a.id, open, changed, !prefix); });
      return add(art, b);
    }
    if (p.closed.length) { add(box, add.apply(null, [el('div', { 'class': 'closed-acts' })].concat(p.closed.map(function (a) { return row(a, false); })))); }
    p.opened.forEach(function (a) { add(box, row(a, true)); });
    return box;
  }
  function toggle(rec, id, open, changed, addressed) {
    var active = records[rec.url] || rec;
    // A control drawn for an act record that has changed since opens nothing: the page draws again instead.
    if (!same(byId(active, id), byId(rec, id))) { changed(); return; }
    var at = addressed ? handles(active).filter(function (a) { return '#f/' + a.handle === location.hash; })[0] : null;
    if (open) { delete active.opened[id]; active.order = active.order.filter(function (x) { return x !== id; }); }
    else { active.opened[id] = true; active.order.push(id); }
    if (at) {
      // The address follows the act it named, under its new handle.
      var now = handles(active).filter(function (a) { return a.id === at.id; })[0];
      history.replaceState(null, '', location.pathname + (now ? '#f/' + now.handle : ''));
    }
    changed();
  }
  function cleanPage() {
    // A closed page has no renderer, search widget, replay, canvas, progress view or external report link.
    document.body.textContent = '';
    Array.prototype.slice.call(root.attributes).forEach(function (a) {
      if (a.name.indexOf('data-pgp-') === 0 || a.name === 'data-finding') { root.removeAttribute(a.name); }
    });
    document.title = 'Results closed: OwnGenome';
  }
  function renderGate() {
    cleanPage();
    // A finding's address names its finding: a closed page shows the run's address instead, the same for any id.
    if (findingPath) { history.replaceState(null, '', new URL('../../', absolute(findingPath)).pathname + location.hash); }
    normalizeHash(current);
    var wrap = el('main', { 'class': 'wrap', id: 'report-gate' });
    add(wrap, el('h1', {}, 'Results about a child'));
    add(wrap, controls(current, function () { if (shut(current)) { renderGate(); } else { transfer(); } }, ''));
    add(document.body, wrap);
    root.setAttribute('data-pgp-status', 'ready');
    root.setAttribute('data-pgp-gated', 'true');
  }
  function proof(rec) {
    var bytes = new TextEncoder().encode(JSON.stringify(rec.acts));
    return crypto.subtle.digest('SHA-256', bytes).then(function (h) { return Array.from(new Uint8Array(h)).map(function (v) { return v.toString(16).padStart(2, '0'); }).join(''); });
  }
  // Replacing the document discards every view closure, image, chart, search cache and in-flight callback. Only
  // which acts are open crosses it, with a digest of their records and a nonce for the next navigation only.
  function transfer() {
    proof(current).then(function (stamp) {
      var ids = current.acts.map(function (f) { return f.id; }), b = new Uint8Array(16);
      crypto.getRandomValues(b);
      var nonce = Array.from(b).map(function (v) { return v.toString(16).padStart(2, '0'); }).join('');
      var target = findingPath || location.pathname;
      sessionStorage.setItem(TRANSFER, JSON.stringify({ path: target, stamp: stamp, nonce: nonce,
        name: window.name.indexOf(NONCE) === 0 ? '' : window.name, order: current.order.map(function (id) { return ids.indexOf(id); }) }));
      window.name = NONCE + nonce;
      location.replace(target); // a map or section address returns to the default view
    });
  }
  function consume(rec) {
    var t = ticket;
    ticket = null;
    if (!t || !rec.acts || !Array.isArray(t.order)) { return Promise.resolve(); }
    return proof(rec).then(function (stamp) {
      if (t.stamp !== stamp) { return; }
      t.order.forEach(function (i) {
        var f = typeof i === 'number' ? rec.acts[i] : null;
        if (f && !rec.opened[f.id]) { rec.opened[f.id] = true; rec.order.push(f.id); }
      });
    });
  }
  // The person and finding pages: the run's state, its family's roles and, for a child, its results, read in the
  // order the views read them (the finding page asks for its results at once, the person page once the state has
  // answered). An adult's results are not read here: the views read them as they would without the gate.
  function preparePerson() {
    var url = absolute(findingPath ? '../../' : './');
    var replaying = root.getAttribute('data-mode') === 'replay';
    var stateRead = sharedJSON(url + (replaying ? 'frames.json' : 'state.json'));
    var results = findingPath ? record(url, undefined, true) : null;
    return stateRead.then(function (raw) {
      var state = replaying ? raw && raw.frames && raw.frames.slice(-1)[0] : raw;
      state = replaying ? state && state.state : state;
      var fam = state && state.plugins && state.plugins.family;
      var entry = fam && fam.kind === 'member' && (fam.families || [])[0];
      return entry ? sharedJSON(absolute(entry.url, url) + 'family.json').then(function (f) {
        var member = f && f.state && (f.state.members || []).filter(function (m) { return m.sample === root.getAttribute('data-run'); })[0];
        // A family that cannot be read, or does not list this person, cannot say they are not a child.
        return !member || !!(member.father || member.mother);
      }) : !state || !state.plugins || !!(fam && fam.unknown);
    }).then(function (child) {
      if (!child && !results) { return { url: url, child: false, acts: [], opened: {}, order: [], perm: [] }; }
      return (results || record(url, child, true)).then(function (rec) { rec.child = child; return rec; });
    }).then(function (rec) {
      current = records[url] = rec;
      return consume(current);
    }).then(function () {
      ownsPage = shut(current);
      if (ownsPage) { renderGate(); }
    }).catch(function () {
      current = { child: true, acts: null, opened: {}, order: [], perm: [] }; ownsPage = true; renderGate();
    });
  }
  // A family's children with a run here; a member who was never sequenced, or has no run, has no report to close.
  function familyRecords(f, base) {
    return Promise.all(((f.state && f.state.members) || []).filter(function (m) {
      return (m.father || m.mother) && f.runs && f.runs[m.sample] && f.runs[m.sample].url;
    }).map(function (m) {
      return record(absolute(f.runs[m.sample].url, base), true).then(function (r) { return [m.sample, r]; });
    })).then(function (entries) { var out = {}; entries.forEach(function (e) { out[e[0]] = e[1]; }); return out; });
  }
  // A run's progress on the home and family pages: how far the run is, never what it found.
  var PROGRESS = ['name', 'kind', 'url', 'state', 'started_utc', 'last_heard_utc', 'finished_utc', 'planned', 'done', 'running', 'failed'];
  function progressOf(r) { var out = {}; PROGRESS.forEach(function (k) { if (r[k] !== undefined) { out[k] = r[k]; } }); return out; }
  function projectFamily(raw, recs) {
    var closed = Object.keys(recs).filter(function (k) { return shut(recs[k]); });
    if (!closed.length) { return raw; }
    function gated(s) { return closed.indexOf(s) >= 0; }
    // An ownership allowlist, not a classifier of genes, prose or individual report fields: a field a producer adds
    // tomorrow stays out until it is listed here.
    var f = { state: { kind: raw.state.kind, name: raw.state.name, pedigree: raw.state.pedigree,
      members: (raw.state.members || []).map(function (m) { return { sample: m.sample, sex: m.sex, father: m.father, mother: m.mother, partner: m.partner, sequenced: m.sequenced }; }),
      children: [], couples: [], plan: [], reports: [] },
      layout: raw.layout ? { width: raw.layout.width, height: raw.layout.height, size: raw.layout.size, lines: copy(raw.layout.lines),
        people: (raw.layout.people || []).map(function (m) { return { sample: m.sample, sex: m.sex, sequenced: m.sequenced, x: m.x, y: m.y }; }) } : null,
      couple_note: raw.couple_note, chrom_lengths: copy(raw.chrom_lengths), gated: closed, runs: {} };
    f.state.children = (raw.state.members || []).filter(function (m) { return m.father || m.mother; }).map(function (m) {
      return gated(m.sample) ? { child: m.sample, gated: true } : copy((raw.state.children || []).filter(function (c) { return c.child === m.sample; })[0]);
    }).filter(Boolean);
    // A row belongs to a gated child when it names that child as one of its members, by exact name.
    function names(v) { return String(v || '').split(/[\s,;]+/); }
    f.state.couples = copy((raw.state.couples || []).filter(function (c) { return !gated(c.a) && !gated(c.b); }));
    // The parents' own items stay: the plan rows that name no gated child, and the couple reports of couples with
    // no gated partner. Any other report of the family folder may hold a child's data and stays out.
    f.state.plan = copy((raw.state.plan || []).filter(function (r) { return !names(r.members).some(gated); }));
    var reports = f.state.couples.map(function (c) { return c.report; }).filter(Boolean);
    f.state.reports = copy((raw.state.reports || []).filter(function (r) { return reports.indexOf(r.url) >= 0; }));
    function free(r) { return !gated(r.sample_a) && !gated(r.sample_b) && !gated(r.a) && !gated(r.b); }
    f.relationships = copy((raw.relationships || []).filter(free));
    f.ibd_pairs = copy((raw.ibd_pairs || []).filter(free));
    f.couples = {};
    Object.keys(raw.couples || {}).forEach(function (key) { if (!key.split(' ').some(gated)) { f.couples[key] = copy(raw.couples[key]); } });
    ['segments', 'trio_status', 'trio_files', 'trio_file_state'].forEach(function (key) {
      if (!raw[key] || typeof raw[key] !== 'object') { return; }
      f[key] = {};
      Object.keys(raw[key]).forEach(function (k) { if (!gated(k)) { f[key][k] = copy(raw[key][k]); } });
    });
    ['de_novo', 'compound_het'].forEach(function (key) { f[key] = copy((raw[key] || []).filter(function (r) { return !gated(r.child); })); });
    Object.keys(raw.runs || {}).forEach(function (k) {
      f.runs[k] = gated(k) ? Object.assign(progressOf(raw.runs[k]), { url: raw.runs[k].url, gated: true }) : copy(raw.runs[k]);
    });
    f.openings = {}; closed.forEach(function (k) { f.openings[k] = (recs[k].acts || []).length; });
    if (raw.side_by_side) {
      var sb = { members: copy(raw.side_by_side.members), rows: [] }, seenRows = Object.create(null);
      // A row is there when an open member has a value in it, never because a gated child has one; the gated
      // child's cell reads "closed". Rows come in a fixed order, so neither a row nor its place depends on the child.
      (raw.side_by_side.rows || []).forEach(function (row) {
        if (!row.values.some(function (v, i) { return !gated(sb.members[i]) && v !== null && v !== undefined && typeof v !== 'object' && ['', '.', '—', 'no summary yet', 'summary not readable'].indexOf(v) < 0; })) { return; }
        var out = { group: row.group, label: row.label,
          values: row.values.map(function (v, i) { return gated(sb.members[i]) ? 'closed' : v; }),
          stale: (row.stale || []).map(function (v, i) { return gated(sb.members[i]) ? false : v; }) };
        var key = JSON.stringify(out);
        if (!seenRows[key]) { seenRows[key] = true; sb.rows.push(out); }
      });
      sb.rows.sort(function (a, b) { var x = JSON.stringify([a.group, a.label]), y = JSON.stringify([b.group, b.label]); return x < y ? -1 : x > y ? 1 : 0; });
      f.side_by_side = sb;
    }
    return f;
  }
  // The home page decides from runs.json alone, as it always has: a run is a child's when a readable family lists
  // one of its parents (or, while a family's members cannot be read, when no readable family lists it), and its
  // summary is closed when runs.json counts or names a result to act on.
  function projectRuns(data) {
    var idx = copy(data), fams = idx.families || [], child = {}, listed = {};
    fams.forEach(function (f) {
      (f.members || []).forEach(function (m) {
        listed[m.sample] = true;
        if ((f.members || []).some(function (x) { return x.sample === m.father || x.sample === m.mother; })) { child[m.sample] = true; }
      });
    });
    var unknown = fams.some(function (f) { return !(f.members || []).length; });
    idx.runs = (idx.runs || []).map(function (r) {
      var n = (r.levels && r.levels.act) || (r.top || []).filter(function (x) { return x && x.level === 'act'; }).length;
      if (!(child[r.name] || (unknown && !listed[r.name])) || !n) { return r; }
      return Object.assign(progressOf(r), { gated: true, results: true, screened: true, attention: null, levels: { act: n },
        top: [{ id: '', level: 'act', title: CLOSED_TEXT }], lead: CLOSED_TEXT, lead_closed: CLOSED_TEXT });
    });
    return idx;
  }
  var familyGates = {}, familyKids = [];
  // The two files the gate projects before a view reads them: family.json on the family page and runs.json on the
  // home page. Every other response reaches its view as it came.
  var familyPage = document.body.classList.contains('page-family');
  function projected(path) { return familyPage ? path.endsWith('/family.json') : path.endsWith('/runs.json'); }
  function loadedJSON(url, data) {
    var u = new URL(url, location.href);
    if (familyPage) {
      familyKids = ((data.state && data.state.children) || []).map(function (c) { return c.child; });
      return familyRecords(data, u.href).then(function (recs) { familyGates = recs; return projectFamily(data, recs); });
    }
    return Promise.resolve(projectRuns(data));
  }
  var ready = Promise.resolve();
  if (document.body.classList.contains('page-person') || document.body.classList.contains('page-finding')) { ready = preparePerson(); }
  // Permission is checked after parsing too: a response started while open may finish after closing.
  function viewResponse(response) {
    var parse = response.json.bind(response);
    response.json = function () { return parse().then(function (data) { return ownsPage ? new Promise(function () {}) : data; }); };
    return response;
  }
  function closedRun(url) { return Object.keys(records).some(function (key) { return url.indexOf(key) === 0 && shut(records[key]); }); }
  // Every view reads through here (and through EventSource below); one that reads a child's closed run gets nothing.
  window.fetch = function (url, options) {
    return ready.then(function () {
      var requestURL = absolute(url && url.url || url);
      if (ownsPage || closedRun(requestURL)) { return new Promise(function () {}); }
      var kept = shared[requestURL];
      delete shared[requestURL];
      var got = kept ? Promise.resolve(new Response(kept.body, { status: kept.status, headers: { 'Content-Type': kept.type } }))
        : originalFetch(url, options);
      return got.then(function (response) {
        if (ownsPage) { return new Promise(function () {}); }
        var path = new URL(requestURL).pathname;
        if (response.ok && current && current.child && path === new URL(current.url + 'results.json').pathname) {
          return response.clone().json().then(function (r) {
            var next = makeRecord(current.url, true, r, current);
            if (!same(next.acts, current.acts)) {
              current = records[current.url] = next;
              // A result to act on that appeared or changed while the report was open closes it at once.
              if (shut(current)) { ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); transfer(); return new Promise(function () {}); }
            }
            return viewResponse(response);
          });
        }
        if (response.ok && projected(path)) {
          return response.json().then(function (data) { return loadedJSON(requestURL, data); }).then(function (data) {
            return viewResponse(new Response(JSON.stringify(data), { status: response.status, headers: response.headers }));
          });
        }
        return viewResponse(response);
      });
    });
  };
  window.addEventListener('hashchange', function () { if (ownsPage) { renderGate(); } });
  // A child's open report never stays in the back-forward cache: it is cleared on the way out, and a page shown
  // again from that cache loads afresh (closed). Other pages come back from the cache as they were.
  window.addEventListener('pagehide', function () {
    if (current && current.child && current.acts && current.acts.length) { ownsPage = true; cleaned = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); }
  });
  window.addEventListener('pageshow', function (e) { if (e.persisted && cleaned) { location.reload(); } });
  var NativeEvents = window.EventSource;
  window.EventSource = function (url, options) {
    var proxy = this, source = null;
    proxy.readyState = 0;
    proxy.close = function () { if (source) { source.close(); } proxy.readyState = 2; };
    ready.then(function () {
      if (ownsPage || closedRun(absolute(url)) || proxy.readyState === 2) { return; }
      source = new NativeEvents(url, options); eventSources.push(source);
      source.onmessage = function (e) { if (!ownsPage && proxy.onmessage) { proxy.onmessage(e); } };
      source.onerror = function (e) { proxy.readyState = source.readyState; if (proxy.onerror) { proxy.onerror(e); } };
      source.onopen = function (e) { proxy.readyState = source.readyState; if (proxy.onopen) { proxy.onopen(e); } };
    });
  };
  window.EventSource.CLOSED = NativeEvents.CLOSED;
  window.PGPUI.ready = ready;
  window.PGPUI.gateOwnsPage = function () { return ownsPage; };
  window.PGPUI.closedText = CLOSED_TEXT;
  // The acts the gate released, for person.js's map switches: every act of a child's open report, else none.
  window.PGPUI.openedActs = function () { return current && !shut(current) ? copy(current.opened) : {}; };
  // Closing one act of a child's open report hands the page back to the gate, which closes it whole.
  window.PGPUI.closeAct = function (id) {
    if (!current || !current.child || !byId(current, id)) { return false; }
    delete current.opened[id];
    current.order = current.order.filter(function (x) { return x !== id; });
    ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); transfer(); return true;
  };
  // What the family page knows of each child's results (for the tests): "DEMO-003:known", "DEMO-005:failed" for a
  // child with no run here.
  window.PGPUI.familyStatus = function () {
    return familyKids.slice().sort().map(function (k) { return k + ':' + (familyGates[k] && familyGates[k].acts ? 'known' : 'failed'); }).join(' ');
  };
  window.PGPUI.familyControls = function (child, redraw) {
    function build() {
      var box = controls(familyGates[child], function () { box.replaceWith(build()); redraw(); }, child + '-');
      box.setAttribute('data-gate-child', child); return box;
    }
    return build();
  };
}());
