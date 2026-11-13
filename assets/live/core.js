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
  // The report gate owns raw payloads. Renderers only see its projection until the final act opens.
  var originalFetch = window.fetch.bind(window);
  var records = {}, current = null, ownsPage = false, eventSources = [];
  var CLOSED_TEXT = "This child's report stays closed until every result to act on in adult life is opened on purpose.";
  var TRANSFER = 'og-report-opening';
  function jsonAt(url) {
    return originalFetch(url, { credentials: 'same-origin', cache: 'no-store' }).then(function (r) {
      if (r.status === 401) { window.location.href = '/login?next=' + encodeURIComponent(location.pathname); }
      return r.ok ? r.json() : null;
    }).catch(function () { return null; });
  }
  function copy(x) { return x === undefined ? null : JSON.parse(JSON.stringify(x)); }
  function absolute(url, base) { return new URL(url, base || location.href).href; }
  function actsOf(r) {
    if (!r || !Array.isArray(r.findings)) { return null; }
    // The producer's act manifest is independent of every other record. Ignore unrelated malformed rows.
    // Level reconciliation supports older documents and refuses to lose an act omitted from a manifest.
    var ids = r.veil ? r.veil.acts : undefined;
    if (ids !== undefined && (!Array.isArray(ids) || ids.some(function (id) { return typeof id !== 'string'; }))) { return null; }
    ids = (ids || []).slice();
    r.findings.forEach(function (f) { if (f && f.level === 'act' && typeof f.id === 'string' && ids.indexOf(f.id) < 0) { ids.push(f.id); } });
    var acts = ids.map(function (id) { return r.findings.filter(function (f) { return f && f.id === id; })[0]; });
    return acts.every(Boolean) ? acts : null;
  }
  function shut(rec) { return rec.child && (!rec.acts || rec.acts.some(function (f) { return !rec.opened[f.id]; })); }
  function record(url, child) {
    url = absolute(url);
    var old = records[url];
    return jsonAt(url + 'results.json').then(function (r) {
      var nextActs = actsOf(r);
      var opened = {};
      (nextActs || []).forEach(function (f) {
        if (old && old.opened[f.id] && (old.acts || []).some(function (prev) { return JSON.stringify(prev) === JSON.stringify(f); })) { opened[f.id] = true; }
      });
      var rec = { url: url, child: child !== undefined ? child : old ? old.child : true, acts: nextActs, opened: opened };
      records[url] = rec;
      return rec;
    });
  }
  function projection(rec) {
    // No data-derived default: even missing/unreadable results have only the fixed explanation.
    return { gated: true, acts: (rec.acts || []).map(function (f) {
      var closed = rec.acts.filter(function (a) { return !rec.opened[a.id]; });
      var opened = rec.acts.filter(function (a) { return rec.opened[a.id]; });
      return rec.opened[f.id] ? { handle: 'opened-' + (opened.indexOf(f) + 1), open: true, finding: copy(f) }
        : { handle: 'closed-' + (closed.indexOf(f) + 1), open: false };
    }) };
  }
  function normalizeHash(p) {
    if (!location.hash && !location.search) { return; }
    var h = location.hash.slice(3);
    var found = /^#f\/(closed|opened)-[1-9][0-9]*$/.test(location.hash) && p.acts.some(function (a) { return a.handle === h; });
    var to = found ? location.hash : p.acts.length ? '#f/' + (p.acts.filter(function (a) { return !a.open; })[0] || p.acts[0]).handle : '';
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
  function controls(rec, changed, only) {
    var p = projection(rec), box = el('section', { 'class': 'report-gate' });
    add(box, el('p', { 'class': 'lede' }, CLOSED_TEXT));
    add(box, el('p', {}, rec.acts ? p.acts.length + ' results to act on' : 'The results cannot be read yet. Reload to try again.'));
    p.acts.forEach(function (a, i) {
      var row = el('article', { id: 'f-' + a.handle, 'data-finding': a.handle, 'class': 'finding' });
      if (a.open && (!only || only === a.finding.id)) { add(row, actRecord(a.finding)); }
      var b = el('button', { type: 'button', 'class': 'btn', 'data-opening': a.handle }, a.open ? 'Close it again' : 'Show the result');
      b.addEventListener('click', function () {
        var active = records[rec.url] || rec, id = rec.acts[i].id;
        var f = (active.acts || []).filter(function (item) { return item.id === id; })[0];
        if (!f || JSON.stringify(f) !== JSON.stringify(rec.acts[i])) { changed(); return; } // stale controls cannot open replacement records
        var at = p.acts.map(function (item) { return '#f/' + item.handle; }).indexOf(location.hash);
        var addressed = at >= 0 ? rec.acts[at].id : null;
        if (a.open) { delete active.opened[f.id]; } else { active.opened[f.id] = true; }
        if (addressed) {
          var next = (active.acts || []).map(function (item) { return item.id; }).indexOf(addressed);
          history.replaceState(null, '', next >= 0 ? '#f/' + projection(active).acts[next].handle : location.pathname);
        }
        changed();
      });
      add(row, b); add(box, row);
    });
    return box;
  }
  function cleanPage() {
    // A gated page has no renderer, search widget, replay, canvas, progress view or external report link.
    document.body.textContent = '';
    Array.prototype.slice.call(root.attributes).forEach(function (a) {
      if (a.name.indexOf('data-pgp-') === 0 || a.name === 'data-finding') { root.removeAttribute(a.name); }
    });
    document.title = 'Results closed: OwnGenome';
    root.removeAttribute('data-pgp-search');
  }
  function renderGate() {
    cleanPage();
    var p = projection(current);
    normalizeHash(p);
    var wrap = el('main', { 'class': 'wrap', id: 'report-gate' });
    add(wrap, el('h1', {}, 'Results about a child'));
    var requested = /\/finding\//.test(location.pathname) ? decodeURIComponent(location.pathname.split('/finding/')[1].split('/')[0]) : null;
    add(wrap, controls(current, function () { if (shut(current)) { renderGate(); } else { transfer(); } }, requested));
    add(document.body, wrap);
    root.setAttribute('data-pgp-status', 'ready');
    root.setAttribute('data-pgp-gated', 'true');
  }
  function proof(rec) {
    var bytes = new TextEncoder().encode(JSON.stringify(rec.acts));
    return crypto.subtle.digest('SHA-256', bytes).then(function (h) { return Array.from(new Uint8Array(h)).map(function (v) { return v.toString(16).padStart(2, '0'); }).join(''); });
  }
  function transfer() {
    // Only opaque booleans and a digest cross the document replacement. The receiver consumes them once.
    // Replacing the document discards every view closure, image, chart, search cache and in-flight callback.
    proof(current).then(function (stamp) {
      sessionStorage.setItem(TRANSFER, JSON.stringify({ path: location.pathname, stamp: stamp, time: Date.now(),
        opened: current.acts.map(function (f) { return !!current.opened[f.id]; }) }));
      location.replace(location.pathname); // requested map/section returns to the default view
    });
  }
  function consume(rec) {
    var text = sessionStorage.getItem(TRANSFER); sessionStorage.removeItem(TRANSFER);
    if (!text) { return Promise.resolve(); }
    return proof(rec).then(function (stamp) {
      try {
        var t = JSON.parse(text);
        if (t.path === location.pathname && t.stamp === stamp && Date.now() - t.time < 5000 && t.opened.length === rec.acts.length) {
          rec.acts.forEach(function (f, i) { if (t.opened[i] === true) { rec.opened[f.id] = true; } });
        }
      } catch (e) { /* An invalid or stale ticket opens nothing. */ }
    });
  }
  function preparePerson() {
    var url = absolute(/\/finding\//.test(location.pathname) ? '../../' : './');
    var replaying = root.getAttribute('data-mode') === 'replay';
    return Promise.all([record(url), jsonAt(url + (replaying ? 'frames.json' : 'state.json'))]).then(function (x) {
      current = x[0];
      var state = replaying ? x[1] && x[1].frames && x[1].frames.slice(-1)[0] : null;
      state = replaying ? state && state.state : x[1];
      var fam = state && state.plugins && state.plugins.family;
      var entry = fam && fam.kind === 'member' && (fam.families || [])[0];
      if (!entry) { current.child = !state || !state.plugins || !!(fam && fam.unknown); return; }
      return jsonAt(absolute(entry.url, url) + 'family.json').then(function (f) {
        var member = f && f.state && (f.state.members || []).filter(function (m) { return m.sample === root.getAttribute('data-run'); })[0];
        current.child = !member || !!(member.father || member.mother);
      });
    }).then(function () { return consume(current); }).then(function () {
      ownsPage = shut(current);
      if (ownsPage) { renderGate(); }
    }).catch(function () {
      current = { child: true, acts: null, opened: {} }; ownsPage = true; renderGate();
    });
  }
  function familyRecords(f, base) {
    return Promise.all(((f.state && f.state.members) || []).filter(function (m) { return m.father || m.mother; }).map(function (m) {
      var run = f.runs && f.runs[m.sample];
      return run && run.url ? record(absolute(run.url, base), true).then(function (r) { return [m.sample, r]; }) : Promise.resolve([m.sample, { child: true, acts: null, opened: {} }]);
    })).then(function (entries) { var out = {}; entries.forEach(function (e) { out[e[0]] = e[1]; }); return out; });
  }
  function projectFamily(raw, recs) {
    var closed = Object.keys(recs).filter(function (k) { return shut(recs[k]); });
    if (!closed.length) { return raw; }
    // An ownership allowlist, not a classifier of genes, prose or individual report fields.
    // Anything a new producer adds stays out until its owner is explicit here.
    var f = { state: { kind: raw.state.kind, name: raw.state.name, pedigree: raw.state.pedigree,
      members: (raw.state.members || []).map(function (m) { return { sample: m.sample, sex: m.sex, father: m.father, mother: m.mother, partner: m.partner, sequenced: m.sequenced }; }), children: [], couples: copy(raw.state.couples || []), plan: [], reports: [] },
      layout: raw.layout ? { width: raw.layout.width, height: raw.layout.height, size: raw.layout.size, lines: copy(raw.layout.lines),
        people: (raw.layout.people || []).map(function (m) { return { sample: m.sample, sex: m.sex, sequenced: m.sequenced, x: m.x, y: m.y }; }) } : null, relationships: copy(raw.relationships || []), ibd_pairs: copy(raw.ibd_pairs || []),
      couples: copy(raw.couples || {}), couple_note: raw.couple_note,
      segments: copy(raw.segments || {}), de_novo: copy(raw.de_novo || []), compound_het: copy(raw.compound_het || []),
      trio_status: copy(raw.trio_status || {}), side_by_side: copy(raw.side_by_side), chrom_lengths: copy(raw.chrom_lengths),
      runs: copy(raw.runs || {}) };
    // Ownership cut: discard child report objects wholesale. Unknown child fields are discarded too.
    f.gated = closed;
    f.openings = {}; closed.forEach(function (k) { f.openings[k] = projection(recs[k]); });
    f.state.children = (raw.state.members || []).filter(function (m) { return m.father || m.mother; }).map(function (m) {
      return closed.indexOf(m.sample) >= 0 ? { child: m.sample, gated: true } : (raw.state.children || []).filter(function (c) { return c.child === m.sample; })[0];
    }).filter(Boolean);
    ['segments', 'trio_status'].forEach(function (key) { closed.forEach(function (k) { if (f[key]) { delete f[key][k]; } }); });
    ['de_novo', 'compound_het'].forEach(function (key) { f[key] = (f[key] || []).filter(function (r) { return closed.indexOf(r.child) < 0; }); });
    f.relationships = (f.relationships || []).filter(function (r) { return !closed.some(function (c) { return JSON.stringify(r).indexOf(c) >= 0; }); });
    f.state.plan = []; f.state.reports = [];
    closed.forEach(function (k) { if (f.runs && f.runs[k]) { f.runs[k] = { url: f.runs[k].url, gated: true }; } });
    f.state.couples = (f.state.couples || []).filter(function (c) { return closed.indexOf(c.a) < 0 && closed.indexOf(c.b) < 0; });
    Object.keys(f.couples || {}).forEach(function (key) { if (key.split(' ').some(function (k) { return closed.indexOf(k) >= 0; })) { delete f.couples[key]; } });
    f.ibd_pairs = (f.ibd_pairs || []).filter(function (r) { return closed.indexOf(r.sample_a) < 0 && closed.indexOf(r.sample_b) < 0; });
    if (f.side_by_side) {
      f.side_by_side = { members: f.side_by_side.members, rows: f.side_by_side.rows };
      var sb = f.side_by_side, seenRows = Object.create(null);
      // Presence comes only from actual open-member values, never an API absence marker.
      sb.rows = sb.rows.filter(function (row) { return row.values.some(function (v, i) { return closed.indexOf(sb.members[i]) < 0 && v !== null && v !== undefined && typeof v !== 'object' && ['', '.', '—', 'no summary yet', 'summary not readable'].indexOf(v) < 0; }); }).map(function (row) {
        row.values = row.values.map(function (v, i) { return closed.indexOf(sb.members[i]) >= 0 ? 'closed' : v; });
        row.stale = (row.stale || []).map(function (v, i) { return closed.indexOf(sb.members[i]) >= 0 ? false : v; });
        return { group: row.group, label: row.label, values: row.values, stale: row.stale };
      }).filter(function (row) {
        var key = JSON.stringify(row);
        if (seenRows[key]) { return false; } seenRows[key] = true; return true;
      }).sort(function (a, b) {
        var x = JSON.stringify([a.group, a.label]), y = JSON.stringify([b.group, b.label]);
        return x < y ? -1 : x > y ? 1 : 0;
      });
    }
    return f;
  }
  var familyGates = {};
  function loadedJSON(url, data) {
    var u = new URL(url, location.href);
    if (u.pathname.endsWith('/family.json') && document.body.classList.contains('page-family')) {
      return familyRecords(data, u.href).then(function (recs) { familyGates = recs; return projectFamily(data, recs); });
    }
    if (u.pathname.endsWith('/runs.json')) {
      var idx = copy(data);
      var families = idx.families || [], listed = {};
      families.forEach(function (f) { (f.members || []).forEach(function (m) { listed[m.sample] = true; }); });
      var unknown = families.some(function (f) { return !(f.members || []).length; });
      var unknownRuns = unknown ? (idx.runs || []).filter(function (r) { return !listed[r.name]; }).map(function (r) {
        return record(absolute(r.url), true).then(function (rec) { var out = {}; out[r.name] = rec; return out; });
      }) : [];
      return Promise.all(unknownRuns.concat(families.map(function (entry) {
        return jsonAt(absolute(entry.url) + 'family.json').then(function (f) {
          if (!f || !f.state) {
            var runs = {};
            (idx.runs || []).forEach(function (r) { runs[r.name] = { url: absolute(r.url) }; });
            f = { state: { members: entry.members || [] }, runs: runs };
          }
          return familyRecords(f, absolute(entry.url));
        });
      }))).then(function (all) {
        all.forEach(function (recs) {
          (idx.runs || []).forEach(function (r, i) {
            if (recs[r.name] && shut(recs[r.name])) {
              idx.runs[i] = { name: r.name, url: r.url, gated: true, state: 'closed', levels: { act: recs[r.name].acts ? recs[r.name].acts.length : 0 },
                top: [{ level: 'act', title: CLOSED_TEXT }], lead: CLOSED_TEXT, lead_closed: CLOSED_TEXT };
            }
          });
        });
        return idx;
      });
    }
    return Promise.resolve(data);
  }
  var ready = Promise.resolve();
  if (document.body.classList.contains('page-person') || document.body.classList.contains('page-finding')) { ready = preparePerson(); }
  // Permission is checked after parsing too: a response started while open may finish after closing.
  function viewResponse(response) {
    var read = response.json.bind(response);
    response.json = function () { return read().then(function (data) { return ownsPage ? new Promise(function () {}) : data; }); };
    return response;
  }
  // All fetch-based readers (including plugins) share this boundary. An unrouted view receives nothing.
  window.fetch = function (url, options) {
    return ready.then(function () {
      var requestURL = absolute(url && url.url || url);
      var closedRun = Object.keys(records).some(function (key) { return requestURL.indexOf(key) === 0 && shut(records[key]); });
      if (ownsPage || closedRun) { return new Promise(function () {}); }
      return originalFetch(url, options).then(function (response) {
        if (ownsPage) { return new Promise(function () {}); }
        var path = new URL(url && url.url || url, location.href).pathname;
        if (response.ok && current && current.child && path === new URL(current.url + 'results.json').pathname) {
          return response.clone().json().then(function (r) {
            var nextActs = actsOf(r);
            if (JSON.stringify(nextActs) !== JSON.stringify(current.acts)) {
              var opened = {};
              (nextActs || []).forEach(function (f) {
                if (current.opened[f.id] && (current.acts || []).some(function (p) { return JSON.stringify(p) === JSON.stringify(f); })) { opened[f.id] = true; }
              });
              current.acts = nextActs; current.opened = opened;
              if (shut(current)) { ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); transfer(); return new Promise(function () {}); }
            }
            return viewResponse(response);
          });
        }
        if (response.ok && (path.endsWith('/family.json') || path.endsWith('/runs.json'))) {
          return response.json().then(function (data) { return loadedJSON(url, data); }).then(function (data) {
            return viewResponse(new Response(JSON.stringify(data), { status: response.status, headers: response.headers }));
          });
        }
        return viewResponse(response);
      });
    });
  };
  window.addEventListener('hashchange', function () { if (ownsPage) { renderGate(); } });
  window.addEventListener('pagehide', function () { if (current && current.child && current.acts && current.acts.length) { ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); } });
  window.addEventListener('pageshow', function (e) { if (e.persisted) { location.reload(); } });
  // live.js also has an event reader. Delay it at the same gate; it never sees child events while closed.
  var NativeEvents = window.EventSource;
  window.EventSource = function (url, options) {
    var proxy = this, source = null;
    proxy.readyState = 0;
    proxy.close = function () { if (source) { source.close(); } proxy.readyState = 2; };
    ready.then(function () {
      var requestURL = absolute(url);
      if (ownsPage || Object.keys(records).some(function (key) { return requestURL.indexOf(key) === 0 && shut(records[key]); }) || proxy.readyState === 2) { return; }
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
  window.PGPUI.openedActs = function () { return current && !shut(current) ? copy(current.opened) : {}; };
  window.PGPUI.closeAct = function (id) {
    if (!current || !current.child || !current.acts || !current.acts.some(function (f) { return f.id === id; })) { return false; }
    delete current.opened[id]; ownsPage = true; eventSources.forEach(function (es) { es.close(); }); cleanPage(); transfer(); return true;
  };
  window.PGPUI.familyStatus = function () { return Object.keys(familyGates).sort().map(function (k) { return k + ':' + (familyGates[k].acts ? 'known' : 'failed'); }).join(' '); };
  window.PGPUI.familyControls = function (child, redraw) {
    function build() {
      var box = controls(familyGates[child], function () { box.replaceWith(build()); redraw(); });
      box.setAttribute('data-gate-child', child); return box;
    }
    return build();
  };

}());
