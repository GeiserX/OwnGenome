// person.js: a person's results on their page (assets/live/index.html), from results.json
// (bin/live_results.py, served by bin/owngenome_live_results.py) and the run's states (live.js).
// The rail holds the report's seven chapter questions, "Look first" (every card worth a look) and
// "What this leaves out"; Progress is live.js's step map. The tab lives in the URL fragment only
// (#q/<chapter>, #s/<section>, #progress), which the browser never sends, so no tab or search reaches
// a server log. A card shows its results once one of its steps is done in the current state, so live
// and replay draw the same way and a replay reveals the results as the run goes.
(function () {
  'use strict';
  if (!window.PGP || !window.PGPUI) { return; }
  var U = window.PGPUI, el = U.el, add = U.add;
  var root = document.documentElement;
  var replay = root.getAttribute('data-mode') === 'replay';
  var results = null, state = null, fetchedFor = null, current = null, chosen = false, drawnKey = null;
  var rail = document.getElementById('rail');
  var panel = document.getElementById('results');
  var progress = document.getElementById('progress');
  var SHORT = { quality: 'Trust', clinical: 'Doctor', carrier: 'Child', medicines: 'Medicines', ancestry: 'Ancestors',
    traits: 'Traits', other: 'Other' };
  var LOOK_EMPTY = 'Nothing on this list. That is information, not a clean bill: each chapter’s caveats and ' +
    '“What this leaves out” still apply.';

  function stepsBy() {
    var by = {};
    ((state && state.steps) || []).forEach(function (s) { by[s.step] = s; });
    return by;
  }
  function doneSet() {
    return ((state && state.steps) || []).filter(function (s) { return s.state === 'done'; })
      .map(function (s) { return s.step; }).join(',');
  }
  function anyDone() { return doneSet() !== ''; }
  // What a card can show at this moment: 'show', 'waiting' (a step still to come), 'notrun' (none planned)
  // or 'empty' (its steps ended and wrote nothing for it).
  function gate(c) {
    var by = stepsBy(), steps = c.steps || [];
    var known = steps.filter(function (s) { return by[s]; });
    if (known.some(function (s) { return by[s].state === 'done'; })) { return c.state === 'missing' ? 'empty' : 'show'; }
    if (known.some(function (s) { return ['waiting', 'running'].indexOf(by[s].state) >= 0; })) { return 'waiting'; }
    if (known.some(function (s) { return by[s].state === 'failed'; })) { return 'failed'; }
    if (!known.length && c.state !== 'missing' && !replay) { return 'show'; }
    return 'notrun';
  }
  function waitingText(c) {
    var by = stepsBy();
    var s = (c.steps || []).filter(function (x) { return by[x] && by[x].state !== 'done'; })[0] || (c.steps || [])[0];
    return 'Waiting for step ' + s + (by[s] ? ' (' + by[s].label + ')' : '') + '.';
  }
  function cards() {
    var out = [];
    ((results && results.chapters) || []).forEach(function (ch) {
      ch.sections.forEach(function (c) { out.push({ ch: ch, c: c }); });
    });
    return out;
  }
  function attentionNow(c) { return c.attention && gate(c) === 'show'; }

  // --- the route --------------------------------------------------------------------------------
  function parse() {
    var h;
    try { h = decodeURIComponent((window.location.hash || '').slice(1)); } catch (e) { return null; }  // a broken fragment
    if (h === 'progress') { return { tab: 'progress' }; }
    if (h.indexOf('q/') === 0) { return { tab: h.slice(2) }; }
    if (h.indexOf('s/') === 0) { return { tab: 'section', id: h.slice(2) }; }
    if (h && cards().some(function (x) { return x.c.id === h; })) { return { tab: 'section', id: h }; }  // an old #clinvar link
    return null;
  }
  function chapterOf(id) {
    var hit = cards().filter(function (x) { return x.c.id === id; })[0];
    return hit ? hit.ch.key : null;
  }
  function defaultTab() { return results && results.chapters && results.chapters.length && anyDone() ? 'look' : 'progress'; }

  function route(focus) {
    var r = parse() || { tab: defaultTab() };
    var tab = r.tab === 'section' ? (chapterOf(r.id) || 'look') : r.tab;
    var valid = ['progress', 'look', 'leaves'].concat(((results && results.chapters) || []).map(function (c) { return c.key; }));
    if (valid.indexOf(tab) < 0) { tab = defaultTab(); }
    show(tab);
    Array.prototype.forEach.call(rail.querySelectorAll('[role=tab]'), function (t) {
      var on = t.getAttribute('data-tab') === tab;
      t.setAttribute('aria-selected', on ? 'true' : 'false');
      t.tabIndex = on ? 0 : -1;
    });
    if (r.tab === 'section') {
      var card = document.getElementById('card-' + r.id);
      if (card) {
        var d = card.querySelector('details.caveat-full');
        if (d && focus !== false) { card.classList.add('target'); }
        card.scrollIntoView({ block: 'start' });
        card.focus({ preventScroll: true });
      }
    } else if (window.location.hash === '#progress') {
      // #progress is also the section's id, so the browser jumps to it: keep the head and the rail in view.
      setTimeout(function () { window.scrollTo(0, 0); }, 0);
    }
  }
  function go(hash) {
    if (window.location.hash === hash) { route(); } else { window.location.hash = hash; }
  }
  window.addEventListener('hashchange', function () { chosen = true; route(); });

  // --- the rail ---------------------------------------------------------------------------------
  function drawRail() {
    rail.textContent = '';
    function tab(key, label, n, hash) {
      var a = el('a', { role: 'tab', href: hash, id: 'tab-' + key, 'data-tab': key, 'aria-controls': key === 'progress' ? 'progress' : 'results' });
      add(a, el('span', { 'class': 'tab-label' }, label));
      if (n) { add(a, el('span', { 'class': 'count', 'aria-label': n + ' worth a look' }, String(n))); }
      add(rail, a);
      return a;
    }
    if (results && results.chapters && results.chapters.length) {
      var all = cards().filter(function (x) { return attentionNow(x.c); }).length;
      tab('look', 'Look first', all, '#q/look');
      results.chapters.forEach(function (ch) {
        var n = ch.sections.filter(attentionNow).length;
        var t = tab(ch.key, SHORT[ch.key] || ch.title, n, '#q/' + ch.key);
        t.setAttribute('title', ch.question);
      });
      tab('leaves', 'Leaves out', 0, '#q/leaves');
    }
    tab('progress', 'Progress', 0, '#progress');
  }
  U.tablist(rail, function (t) { go(t.getAttribute('href')); });

  // --- a card -----------------------------------------------------------------------------------
  function badge(text, cls, symbol) {
    var b = el('span', { 'class': 'chip ' + cls });
    if (symbol) { add(b, el('span', { 'class': 'sym', 'aria-hidden': 'true' }, symbol)); }
    add(b, text);
    return b;
  }
  function recheckBadges(c) {
    var r = state && state.plugins && state.plugins.recheck;
    var out = [];
    (c.steps || []).forEach(function (s) {
      ((r && r.steps && r.steps[s]) || []).forEach(function (t) { out.push(badge(t, 'chip-recheck', '↻')); });
    });
    return out;
  }
  function card(c, ch) {
    var g = gate(c);
    var cls = 'card result g-' + g + (attentionNow(c) ? ' attn' : '');
    var art = el('article', { 'class': cls, id: 'card-' + c.id, 'data-section': c.id, 'data-step': (c.steps || []).join(' '),
      'data-state': c.state, tabindex: '-1', 'aria-labelledby': 'h-' + c.id });
    var head = el('header', { 'class': 'card-head' });
    var h = el('h3', { id: 'h-' + c.id });
    add(h, el('a', { href: '#s/' + c.id, 'class': 'self' }, c.title));
    add(head, h);
    var chips = el('div', { 'class': 'chips' });
    if (g === 'show') {
      add(chips, attentionNow(c) ? badge('worth a look', 'chip-attn', '▲') : badge('for information', 'chip-info', '○'));
    }
    recheckBadges(c).forEach(function (b) { add(chips, b); });
    add(head, chips);
    add(art, head);
    if (g === 'waiting') {
      add(art, el('p', { 'class': 'state-text' }, waitingText(c)));
    } else if (g === 'notrun') {
      add(art, el('p', { 'class': 'state-text' }, 'Not run in this run (step ' + (c.steps || []).join(', ') + ').'));
    } else if (g === 'failed') {
      add(art, el('p', { 'class': 'state-text' }, 'Step ' + (c.steps || []).join(', ') + ' failed: the Progress tab says why.'));
    } else if (g === 'empty') {
      add(art, el('p', { 'class': 'state-text' }, 'Its step ended without a result for this section.'));
    } else {
      (c.notes || []).forEach(function (n) { add(art, el('p', { 'class': 'note-stale' }, n)); });
      if (c.headline) { add(art, el('p', { 'class': 'headline' }, c.headline)); }
      add(art, U.facts(c.facts));
      Object.keys((results && results.charts) || {}).forEach(function (chart) {
        if (results.charts[chart] === c.id) { add(art, el('div', { 'class': 'chart-slot', 'data-chart-slot': chart })); }
      });
      if (c.image) { add(art, el('img', { 'class': 'plot', src: c.image.src, alt: c.image.alt })); }
      var tb = U.table(c.table);
      if (tb) { add(art, tb); } else if (c.empty_text) { add(art, el('p', { 'class': 'empty' }, c.empty_text)); }
      if (c.lines && c.lines.length) {
        var ul = el('ul', { 'class': 'lines' });
        c.lines.forEach(function (l) { add(ul, el('li', {}, l)); });
        add(art, ul);
      }
      if (c.hidden) {
        var d = el('details', { 'class': 'behind' });
        add(d, el('summary', {}, c.hidden.summary), el('p', {}, c.hidden.text));
        add(art, d);
      }
      if (c.explore && c.explore.length) {
        var p = el('p', { 'class': 'explore' }, 'Explore the branch (opens a website; the address carries the haplogroup): ');
        c.explore.forEach(function (x, i) {
          if (i) { add(p, ', '); }
          add(p, el('a', { href: x.url, rel: 'noreferrer noopener', target: '_blank' }, x.name));
        });
        add(art, p);
      }
    }
    add(art, U.caveat(c.caveat));
    var foot = el('p', { 'class': 'card-links' });
    add(foot, el('a', { href: 'report.html#' + encodeURIComponent(c.anchor), 'class': 'report' }, 'Open in the report'));
    if (g === 'show' && c.viewer && c.viewer.length) {
      add(foot, el('a', { href: 'viewer/#' + c.viewer[0], 'class': 'viewer' }, 'See the reads at ' + c.viewer[0]));
    }
    (g === 'show' ? c.links || [] : []).forEach(function (f) {
      add(foot, el('a', { href: f.split('/').map(encodeURIComponent).join('/') }, 'Open ' + f));
    });
    add(art, foot);
    return art;
  }

  // --- the views --------------------------------------------------------------------------------
  function heading(text, sub) {
    var h = el('header', { 'class': 'view-head' });
    add(h, el('h2', { 'class': 'question' }, text));
    if (sub) { add(h, el('p', { 'class': 'dim' }, sub)); }
    return h;
  }
  function lookFirst() {
    var box = el('div', { id: 'look-first', 'class': 'view' });
    add(box, heading('Look at these first', 'The cards the report itself flags, in chapter order. A card not listed is for information; its caveats still apply.'));
    var items = cards().filter(function (x) { return attentionNow(x.c); });
    if (!items.length) {
      add(box, el('p', { 'class': 'empty look-empty' }, LOOK_EMPTY));
    } else {
      var ol = el('ol', { 'class': 'look-list' });
      items.forEach(function (x) {
        var li = el('li', {});
        var a = el('a', { href: '#s/' + x.c.id });
        add(a, el('span', { 'class': 'look-ch' }, SHORT[x.ch.key] || x.ch.title), el('span', { 'class': 'look-title' }, x.c.title));
        add(li, a);
        if (x.c.headline) { add(li, el('p', { 'class': 'look-head' }, x.c.headline)); }
        add(ol, li);
      });
      add(box, ol);
    }
    var waiting = cards().filter(function (x) { return gate(x.c) === 'waiting'; }).length;
    if (waiting) { add(box, el('p', { 'class': 'dim' }, (waiting === 1 ? '1 card is' : waiting + ' cards are') + ' still waiting for its step.')); }
    return box;
  }
  function chapter(ch) {
    var box = el('div', { 'class': 'view chapter', 'data-chapter': ch.key });
    add(box, heading(ch.question, ch.title));
    var grid = el('div', { 'class': 'cards' });
    ch.sections.forEach(function (c) { add(grid, card(c, ch)); });
    add(box, grid);
    return box;
  }
  function leaves() {
    var box = el('div', { 'class': 'view' });
    add(box, heading('What this leaves out', 'What no step looks at, and the steps that did not run.'));
    function list(title, items) {
      var s = el('section', { 'class': 'card' });
      add(s, el('h3', {}, title));
      var ul = el('ul', { 'class': 'lines' });
      (items.length ? items : ['None.']).forEach(function (t) { add(ul, el('li', {}, t)); });
      add(s, ul);
      add(box, s);
    }
    list('No step looks at these', (results && results.not_assessed) || []);
    list('Steps that did not run', (results && results.not_run) || []);
    return box;
  }
  function show(tab) {
    current = tab;
    var isProgress = tab === 'progress';
    progress.hidden = !isProgress;
    panel.hidden = isProgress;
    root.setAttribute('data-pgp-tab', tab);
    if (isProgress) { return; }
    panel.textContent = '';
    if (results && results.note && !(results.chapters || []).length) {
      add(panel, el('p', { 'class': 'empty' }, results.note));
      return;
    }
    if (tab === 'look') { add(panel, lookFirst()); } else if (tab === 'leaves') { add(panel, leaves()); } else {
      var ch = ((results && results.chapters) || []).filter(function (c) { return c.key === tab; })[0];
      if (ch) { add(panel, chapter(ch)); }
    }
    root.setAttribute('data-pgp-results', String(panel.querySelectorAll('article.result').length));
    document.dispatchEvent(new CustomEvent('pgp:cards'));
  }

  // --- the header -------------------------------------------------------------------------------
  function header() {
    var p = (results && results.person) || {};
    var slot = document.getElementById('avatar-slot');
    var planned = ((state && state.steps) || []).filter(function (s) { return s.state !== 'not planned'; });
    var done = planned.filter(function (s) { return s.state === 'done'; }).length;
    slot.textContent = '';
    add(slot, U.avatar(p.inferred_sex || (state && state.run && state.run.declared_sex), true, planned.length ? done / planned.length : 0));
    var line = document.getElementById('person-line');
    line.textContent = '';
    var bits = [];
    if (p.inferred_sex) { bits.push('Sex from the reads: ' + p.inferred_sex + (p.declared_sex ? ', declared ' + p.declared_sex : '')); }
    if (bits.length) { add(line, bits.join('. ') + '. '); }
    var fam = state && state.plugins && state.plugins.family;
    var crumbs = document.getElementById('crumbs');
    var name = document.getElementById('run-name');
    crumbs.textContent = '';
    ((fam && fam.kind === 'member' && fam.families) || []).forEach(function (f, i) {
      if (i === 0) { add(line, 'In family '); } else { add(line, ', '); }
      add(line, el('a', { href: f.url }, f.name));
      add(crumbs, el('a', { href: f.url }, f.name), el('span', { 'class': 'sep', 'aria-hidden': 'true' }, '/'));
    });
    add(crumbs, name);
    // The strip: one cell per planned step, filled when done, striped while it runs.
    var strip = document.getElementById('strip');
    strip.textContent = '';
    planned.forEach(function (s) {
      add(strip, el('span', { 'class': 'cell st-' + s.state.replace(/\s+/g, '-'), title: s.step + ' ' + s.label + ': ' + s.state }));
    });
    var txt = document.getElementById('strip-text');
    var running = planned.filter(function (s) { return s.state === 'running'; }).length;
    var failed = planned.filter(function (s) { return s.state === 'failed'; }).length;
    var run = (state && state.run) || {};
    var t = done + ' of ' + planned.length + ' steps done';
    if (run.state === 'finished') { t = 'Finished' + (run.finished_utc ? ' ' + U.day(run.finished_utc) : '') + ': ' + t; }
    else if (running) { t = 'Running: ' + t + ', ' + running + ' running now'; }
    if (failed) { t += ', ' + failed + ' failed'; }
    var rc = state && state.plugins && state.plugins.recheck;
    if (rc && rc.checked_utc) { t += '. Databases re-checked ' + U.day(rc.checked_utc); }
    txt.textContent = t + '.';
  }

  // --- data ---------------------------------------------------------------------------------------
  function load() {
    var key = replay ? 'replay' : doneSet() + '|' + ((state && state.results_stamp) || '');
    if (key === fetchedFor) { return; }
    fetchedFor = key;
    U.getJSON('results.json').then(function (r) {
      if (key !== fetchedFor) { return; }
      if (!(r && (r.chapters || r.note))) { fetchedFor = null; return; }  // tried again on the next state
      results = r;
      redraw(true);
    }).catch(function () { if (key === fetchedFor) { fetchedFor = null; } });
  }
  function redraw(force) {
    var k = JSON.stringify([doneSet(), state && state.run && state.run.state, !!results,
      state && state.plugins && state.plugins.recheck && state.plugins.recheck.mtime]);
    if (!force && k === drawnKey) { return; }
    drawnKey = k;
    header();
    var keep = document.activeElement && rail.contains(document.activeElement) ? document.activeElement.getAttribute('data-tab') : null;
    drawRail();
    if (!chosen && !parse()) { show(defaultTab()); route(false); } else { route(false); }
    if (keep) { var t = rail.querySelector('[data-tab="' + keep + '"]'); if (t) { t.focus(); } }
  }

  window.PGP.register({
    name: 'person',
    onState: function (s) {
      state = s;
      if (!replay || !results) { load(); }
      redraw(false);
    }
  });
}());
