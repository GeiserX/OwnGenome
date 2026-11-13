// person.js: a person's genome as a report, the Journal design (docs/design/journal/), on their page
// (assets/live/index.html), from results.json (bin/live_results.py, served by
// bin/owngenome_live_results.py) and the run's states (live.js). The page reads top to bottom: the lead
// paragraph built from the findings, the counts by level, the genome map of the findings (map.js), then ten
// numbered sections: 1 To act on, 2 Medicines, 3 For your family, 4 Worth a look, 5 Can I trust this?,
// 6 Origins, 7 Traits and scores, 8 Assessed, no finding, 9 What it cannot see, 10 Sources and files.
// A run still going is the same report being written: each section carries its stamp (written, being
// written, waiting), a short band under the header counts the steps, and the step timeline closes the
// page, folded. Every card of results.json is on the page under
// its section; a card shows its results once one of its steps is done in the current state, so live
// and replay draw the same way. Links inside the page live in the fragment (#s/<card>, #s3, #f/<id>,
// #progress), which the browser never sends. Silence is never an all-clear: a section with no finding
// says what its steps looked at and what they do not see.
(function () {
  'use strict';
  if (!window.PGP || !window.PGPUI) { return; }
  var U = window.PGPUI, el = U.el, add = U.add;
  var root = document.documentElement;
  var replay = root.getAttribute('data-mode') === 'replay';
  var results = null, state = null, family = null, familyFor = null, fetchedFor = null, drawnKey = null, folded = false;
  // 'none': in no family; 'pending' or 'failed': in one whose roles are not known (yet); 'known'.
  var familyStatus = 'none';
  var shown = {};  // veiled findings the reader opened, by id
  var genomeMap = null, resultsKey = null;  // Figure 1 (map.js), and the results it was built for
  var panel = document.getElementById('results');
  var rail = document.getElementById('rail');
  var NOT_ALL_CLEAR = 'No finding in what these steps assessed.';

  function stepsBy() {
    var by = {};
    ((state && state.steps) || []).forEach(function (s) { by[s.step] = s; });
    return by;
  }
  function doneSet() {
    return ((state && state.steps) || []).filter(function (s) { return s.state === 'done'; })
      .map(function (s) { return s.step; }).join(',');
  }
  function running() { return !!(state && state.run && state.run.state === 'running'); }
  // What a card can show at this moment: 'show', 'waiting' (a step still to come), 'running' (one of its
  // steps runs now), 'failed', 'notrun' (none planned) or 'empty' (its steps ended and wrote nothing).
  function gate(c) {
    var by = stepsBy(), steps = c.steps || [];
    var known = steps.filter(function (s) { return by[s]; });
    if (known.some(function (s) { return by[s].state === 'done'; })) { return c.state === 'missing' ? 'empty' : 'show'; }
    if (known.some(function (s) { return by[s].state === 'running'; })) { return 'running'; }
    if (known.some(function (s) { return by[s].state === 'waiting'; })) { return 'waiting'; }
    if (known.some(function (s) { return by[s].state === 'failed'; })) { return 'failed'; }
    if (!known.length && c.state !== 'missing' && !replay) { return 'show'; }
    return 'notrun';
  }
  function pending(g) { return g === 'waiting' || g === 'running'; }
  function stepWords(c) {
    var by = stepsBy();
    var s = (c.steps || []).filter(function (x) { return by[x] && by[x].state !== 'done'; })[0] || (c.steps || [])[0];
    return s ? (by[s] ? by[s].label : 'step ' + s) : 'its step';
  }
  function cards() {
    var out = [];
    ((results && results.chapters) || []).forEach(function (ch) { ch.sections.forEach(function (c) { out.push(c); }); });
    return out;
  }
  function cardById(id) { return cards().filter(function (c) { return c.id === id; })[0]; }
  function findingById(id) { return ((results && results.findings) || []).filter(function (f) { return f.id === id; })[0]; }
  function findingsNow() {
    return ((results && results.findings) || []).filter(function (f) { var c = cardById(f.card); return c && gate(c) === 'show'; });
  }
  function sectionsList() { return (results && results.journal) || []; }
  function attentionNow(c) { return c.attention && gate(c) === 'show'; }
  function personName() { return (results && results.person && results.person.sample) || root.getAttribute('data-run') || ''; }

  // --- the family this person is in: family.json of the first family the family plugin names -------
  function familyUrl() {
    var fam = state && state.plugins && state.plugins.family;
    var f = fam && fam.kind === 'member' && (fam.families || [])[0];
    return f ? f : null;
  }
  // The family plugin cannot say: a state with no plugins (unavailable or cut short), or a family plugin
  // that found a family report it could not read (plugins.family.unknown). Not a standalone adult then.
  function familyUnknown() {
    var fam = state && state.plugins ? state.plugins.family : undefined;
    return !state || !state.plugins || !!(fam && fam.unknown);
  }
  function loadFamily() {
    var f = familyUrl();
    var key = f ? f.url : (familyUnknown() ? '?' : null);
    if (key === familyFor) { return; }
    familyFor = key;
    if (!f) { family = null; familyStatus = key ? 'failed' : 'none'; return; }
    familyStatus = 'pending';
    U.getJSON(f.url + 'family.json').then(function (d) {
      if (familyFor !== key) { return; }
      family = d && d.state ? { name: f.name, url: f.url, data: d } : null;
      familyStatus = family ? 'known' : 'failed';
      redraw(true);
    }).catch(function () {
      if (familyFor !== key) { return; }  // a late answer for a family this page no longer shows
      family = null;
      familyStatus = 'failed';
      redraw(true);
    });
  }
  function me() {
    var ms = (family && family.data.state && family.data.state.members) || [];
    return ms.filter(function (m) { return m.sample === personName(); })[0] || null;
  }
  // A child of the family: a member whose father or mother is listed. Results to act on about a child
  // stay closed until someone opens them on purpose.
  function isChild() { var m = me(); return !!(m && (m.father || m.mother)); }
  // Results to act on stay closed for a child, and for anyone in a family whose roles are not known yet:
  // the page cannot tell a child apart until family.json has answered.
  // A family that answers without this person cannot say they are not a child either.
  function closedPage() { return isChild() || familyStatus === 'pending' || familyStatus === 'failed' || (familyStatus === 'known' && !me()); }
  function memberUrl(s) {
    var r = family && family.data.runs && family.data.runs[s];
    return r && r.url ? r.url : null;
  }
  function role() {
    var m = me();
    if (!m) { return null; }
    var ms = family.data.state.members;
    var kids = ms.filter(function (x) { return x.father === m.sample || x.mother === m.sample; });
    var bits = [];
    if (m.father || m.mother) {
      bits.push(({ male: 'son', female: 'daughter' }[m.sex] || 'child') + ' of ' + [m.father, m.mother].filter(Boolean).join(' and '));
    }
    if (kids.length) { bits.push(({ male: 'father', female: 'mother' }[m.sex] || 'parent') + ' of ' + kids.map(function (k) { return k.sample; }).join(', ')); }
    if (m.partner) { bits.push('partner of ' + m.partner); }
    return bits.join(', ');
  }

  // --- widgets: a card's numbers drawn ---------------------------------------------------------------
  // HTML placed in percent of a track (custom properties, never a style attribute), so the text stays
  // at its size on a phone. Every mark is focusable and says itself in words; colour is never alone.
  function num(v) { var n = parseFloat(v); return isFinite(n) ? n : null; }
  function pct(v, max) { return Math.max(0, Math.min(100, 100 * v / max)); }
  function nice(v) { var step = v > 200 ? 50 : 10; return Math.max(step, Math.ceil(v / step) * step); }
  function place(e, from, to, max) {
    e.style.setProperty('--a', pct(from, max).toFixed(2) + '%');
    if (to !== undefined) { e.style.setProperty('--w', Math.max(0.6, pct(to, max) - pct(from, max)).toFixed(2) + '%'); }
    return e;
  }
  function mark(cls, v, max, label) {
    var m = el('span', { 'class': 'w-tick ' + (cls || ''), tabindex: '0', role: 'img', 'aria-label': label, title: label });
    return place(m, v, undefined, max);
  }
  function trackOf(max, parts) {
    var t = el('div', { 'class': 'w-track' });
    parts.forEach(function (p) { if (p) { add(t, p); } });
    return t;
  }
  function band(cls, from, to, max, label) {
    return place(el('span', { 'class': 'w-band ' + cls, title: label || '' }), from, to, max);
  }
  function whisker(lo, hi, max, label) {
    return place(el('span', { 'class': 'w-ci', title: label || '' }), lo, hi, max);
  }
  function wrow(head, track, legend, cls) {
    var r = el('div', { 'class': 'w-row ' + (cls || '') });
    add(r, head, track);
    if (legend) { add(r, el('p', { 'class': 'w-legend' }, legend)); }
    return r;
  }
  function whead(name, value, extra) {
    var h = el('p', { 'class': 'w-head' });
    add(h, el('span', { 'class': /^[A-Z0-9-]+( \(|$)/.test(name) ? 'w-name' : 'w-label' }, name), ' ', el('span', { 'class': 'w-val' }, value));
    if (extra) { add(h, ' ', el('span', { 'class': 'w-extra' }, extra)); }
    return h;
  }

  // Repeat expansions: each allele against Stranger's thresholds, its 95% interval as a whisker.
  function strGauge(w) {
    var box = el('div', { 'class': 'widget w-str', 'data-widget': 'str' });
    (w.loci || []).forEach(function (l) {
      var counts = String(l.repeat_count || '').split('/');
      var sizing = String(l.allele_sizing || '').split('/');
      var cis = String(l.ci || '').split('/');
      var nmax = num(l.normal_max), pmin = num(l.pathologic_min);
      var top = 0;
      counts.forEach(function (c, i) {
        var ci = (cis[i] || '').split('-');
        top = Math.max(top, num(c) || 0, num(ci[1]) || 0);
      });
      var max = nice(Math.max(top * 1.15, (pmin || 0) * 1.3, (nmax || 0) * 1.6, 10));
      var parts = [];
      if (nmax !== null) { parts.push(band('w-normal', 0, nmax, max, 'Stranger’s "normal": up to ' + nmax)); }
      if (pmin !== null) { parts.push(band('w-patho', pmin, max, max, 'pathologic from ' + pmin)); }
      counts.forEach(function (c, i) {
        var v = num(c);
        if (v === null) { return; }
        var ci = (cis[i] || '').split('-');
        var lo = num(ci[0]), hi = num(ci[1]);
        if (lo !== null && hi !== null && hi > lo) { parts.push(whisker(lo, hi, max, '95% interval ' + lo + ' to ' + hi)); }
        var unsized = sizing[i] && ['sized', 'unknown'].indexOf(sizing[i]) < 0;
        parts.push(mark(unsized ? 'w-low' : '', v, max, l.locus + ' allele ' + (i + 1) + ': ' + c + ' repeats' +
          (lo !== null && hi !== null ? ', 95% interval ' + lo + ' to ' + hi : '') + (unsized ? ', not sized from short reads' : '')));
      });
      var st = l.status ? 'Stranger: "' + l.status + '"' : 'no Stranger status';
      var legend = (nmax !== null || pmin !== null)
        ? [nmax !== null ? '"normal" up to ' + nmax : '', pmin !== null ? 'pathologic from ' + pmin : '', 'scale 0 to ' + max]
          .filter(Boolean).join(' · ')
        : 'No Stranger thresholds for this locus. Scale 0 to ' + max + '.';
      var flagged = l.status && l.status.split(',').some(function (x) { return x && x !== 'normal'; });
      add(box, wrow(whead(l.locus, (l.repeat_count || '.') + ' repeats', st), trackOf(max, parts), legend, flagged ? 'w-attn' : ''));
    });
    add(box, el('p', { 'class': 'w-key dim' }, 'Tick: an allele’s repeat count; line: its 95% interval; hollow: not sized ' +
      'from short reads. Shaded: Stranger’s "normal" band; hatched: its pathologic range.'));
    return box;
  }

  // A count with its 95% interval on a scale from 0 (Lp(a)'s KIV-2 repeats).
  function ciBar(w) {
    var box = el('div', { 'class': 'widget w-cibar', 'data-widget': 'ci' });
    var max = w.max || nice(w.high * 1.25);
    var label = w.label + ': ' + w.value + ', 95% interval ' + w.low + ' to ' + w.high;
    add(box, wrow(whead(w.label, w.value + ' ' + (w.unit || ''), '95% interval ' + w.low + ' to ' + w.high +
      (w['class'] ? ' · class "' + w['class'] + '"' : '')),
      trackOf(max, [whisker(w.low, w.high, max, '95% interval'), mark('', w.value, max, label)]),
      'Scale 0 to ' + max + ' ' + (w.unit || '') + '. The line is the 95% interval of the estimate.'));
    return box;
  }

  // CHIP: each call's VAF against the smallest VAF this depth can see in its own gene. The genes differ:
  // below the smallest limit no gene sees a clone, between the smallest and the largest some genes do.
  function vafChart(w) {
    var box = el('div', { 'class': 'widget w-vaf', 'data-widget': 'vaf' });
    var max = w.max || 0.5;
    var limits = (w.genes || []).map(function (g) { return g.limit; }).filter(function (x) { return typeof x === 'number'; });
    var lo = limits.length ? Math.min.apply(null, limits) : null, hi = limits.length ? Math.max.apply(null, limits) : null;
    var parts = [];
    if (lo !== null) { parts.push(band('w-blind', 0, lo, max, 'below the detection limit of every gene at this depth')); }
    if (hi !== null && hi > lo) { parts.push(band('w-part', lo, hi, max, 'below the detection limit of some genes at this depth')); }
    // A call's own gene limit, in words: a call below it is a clone that depth sees only now and then.
    var own = function (c) {
      return typeof c.limit === 'number' ? (c.vaf < c.limit ? 'below' : 'above') + ' the limit of ' + c.gene + ' at this depth (' + c.limit + ')'
        : 'no detection limit for ' + c.gene;
    };
    (w.calls || []).forEach(function (c) {
      parts.push(mark('w-call' + (typeof c.limit === 'number' && c.vaf < c.limit ? ' w-low' : ''), c.vaf, max,
        c.gene + ' ' + (c.change || '') + ': VAF ' + c.vaf + ', ' + own(c)));
    });
    var head = whead('Variant allele fraction', (w.calls || []).length ? (w.calls.length + ' call' + (w.calls.length > 1 ? 's' : ''))
      : 'no call above the limit');
    var legend = (lo !== null ? 'Hatched: below the smallest VAF this depth can see in any of the ' + limits.length + ' genes (' + lo +
      '); a clone that small is not seen. ' + (hi > lo ? 'Shaded: up to ' + hi + ', seen in some genes and not in others (the table below has each gene\'s limit). ' : '') : '') +
      'Scale VAF 0 to ' + max + '.';
    add(box, wrow(head, trackOf(max, parts), legend));
    if ((w.calls || []).length) {
      var ul = el('ul', { 'class': 'w-calls' });
      w.calls.forEach(function (c) { add(ul, el('li', {}, c.gene + ' ' + (c.change || '') + ' at VAF ' + c.vaf + ', ' + own(c))); });
      add(box, ul);
    }
    if ((w.genes || []).length) {
      var d = el('details', { 'class': 'w-fallback' });
      add(d, el('summary', {}, 'Detection limit per gene'));
      add(d, U.table({ columns: ['Gene', 'Mean depth', 'Smallest VAF seen', 'Status'],
        rows: w.genes.map(function (g) { return [g.gene, g.depth === null ? '.' : String(g.depth), g.limit === null ? '.' : String(g.limit), g.status || '.']; }) }, 40));
      add(box, d);
    }
    return box;
  }

  // Copy numbers as dots, one per copy; a call that is not "called" is hollow and says so.
  function copyDots(w) {
    var box = el('div', { 'class': 'widget w-copies', 'data-widget': 'copies' });
    var ul = el('ul', { 'class': 'w-copy-list' });
    (w.genes || []).forEach(function (g) {
      var low = g.status !== 'called';
      var li = el('li', { 'class': low ? 'w-low' : '' });
      add(li, el('span', { 'class': 'w-name' }, g.gene));
      var dots = el('span', { 'class': 'w-dots', role: 'img', tabindex: '0',
        'aria-label': g.gene + ': ' + (g.copies === null ? 'no copy number' : g.copies + ' copies') + (low ? ', ' + g.status : '') });
      for (var i = 0; i < Math.min(g.copies || 0, 8); i++) { add(dots, el('span', { 'class': 'w-dot' })); }
      if ((g.copies || 0) > 8) { add(dots, '+' + (g.copies - 8)); }
      add(li, dots);
      add(li, el('span', { 'class': 'w-val' }, g.copies === null ? 'not called' : g.copies + (g.copies === 1 ? ' copy' : ' copies')));
      if (low) { add(li, el('span', { 'class': 'w-extra' }, 'low confidence: ' + g.status)); }
      add(ul, li);
    });
    add(box, ul);
    return box;
  }

  // The Y haplogroup's path from the root of the tree to its placement.
  function yPath(w) {
    var box = el('div', { 'class': 'widget w-ypath', 'data-widget': 'ypath' });
    add(box, el('p', { 'class': 'w-head' }, 'From the root to ' + (w.tip || 'the placement') + (w.tree ? ' on ' + w.tree : '')));
    var ol = el('ol', { 'class': 'w-path' });
    (w.path || []).forEach(function (p, i, all) {
      var seen = (p.derived || 0) + (p.ancestral || 0) + (p.no_call || 0);
      var li = el('li', { 'class': (i === all.length - 1 ? 'w-tip' : '') + (p.ancestral ? ' w-low' : ''), tabindex: '0',
        'aria-label': 'Level ' + p.level + ', ' + p.branch + ': ' + p.derived + ' of ' + seen + ' markers derived' +
          (p.ancestral ? ', ' + p.ancestral + ' ancestral' : '') });
      add(li, el('span', { 'class': 'w-branch' }, p.branch));
      add(li, el('span', { 'class': 'w-count' }, (p.derived === null ? '.' : p.derived) + ' of ' + seen + ' derived' +
        (p.ancestral ? ' · ' + p.ancestral + ' ancestral' : '')));
      add(ol, li);
    });
    add(box, ol);
    return box;
  }
  var WIDGETS = { str: strGauge, ci: ciBar, vaf: vafChart, copies: copyDots, ypath: yPath };
  function widget(w) {
    var fn = w && WIDGETS[w.kind];
    if (!fn) { return null; }
    try { return fn(w); } catch (e) { return el('p', { 'class': 'dim' }, 'This drawing could not be made; the table below has the numbers.'); }
  }

  // --- a card: one section of the report, with its numbers, its rows, its figure and its caveat -------
  function recheckWords(c) {
    var r = state && state.plugins && state.plugins.recheck;
    var out = [];
    (c.steps || []).forEach(function (s) {
      ((r && r.steps && r.steps[s]) || []).forEach(function (t) { if (out.indexOf(t) < 0) { out.push(t); } });
    });
    return out;
  }
  function chip(text, cls, symbol) {
    var b = el('span', { 'class': 'chip ' + cls });
    if (symbol) { add(b, el('span', { 'class': 'sym', 'aria-hidden': 'true' }, symbol)); }
    add(b, text);
    return b;
  }
  // The findings of a card that a child's page keeps closed: while one is closed the card's own rows,
  // numbers and drawing stay closed too, or its table would name the gene and variant under the veil.
  // A place outside its result's card (a ClinVar row, a card's line, fact or headline, a viewer link, a chart's mark,
  // a search hit) carries a tag (results.veil, bin/live_results.py veil()): the list of every result to act on it
  // names. It stays closed while ANY of them is, so a place naming two closed results never comes back with the
  // first one opened. A tag naming no finding of these results stays closed. PANEL (a fixed panel's line, or prose
  // that cannot be shown to name no gene) names no result: it stays closed while any result to act on is,
  // whichever, so what is missing never tells the gene. The ids of a tag never reach the page.
  var PANEL = '*panel*';
  function tagClosed(id) {
    if (Array.isArray(id)) { return id.some(tagClosed); }
    if (id === PANEL) { return closedPage() && closedActs().length > 0; }
    var f = findingById(id);
    return closedPage() && (!f || veiled(f));
  }
  // Explicit null is a structured place with no matching act. Missing metadata is unclassified: PANEL.
  function placeClosed(tag) { return tag === undefined ? tagClosed(PANEL) : !!tag && tagClosed(tag); }
  function veilTags(kind, key) { return ((results && results.veil && results.veil[kind]) || {})[key] || []; }
  // A card's headline is prose: tagged PANEL unless the summary key it comes from prints no gene (veil.headline).
  function headlineShut(c) { var t = ((results && results.veil && results.veil.headline) || {})[c.id]; return placeClosed(t); }
  // The charts (plugins/charts.js) ask for each mark. A closed page with no results.json yet cannot tell which
  // marks those are, so it leaves every mark out until it can.
  window.PGP.markClosed = function (chart, mark) {
    if (!closedPage()) { return false; }
    if (!results) { return true; }
    var id = ((results.veil && results.veil.marks && results.veil.marks[chart]) || {})[mark];
    return placeClosed(id);
  };
  // Chart prose and dimensions have no structured result tag. The chart reader asks the same PANEL rule.
  function chartPanelClosed() { return closedPage() && (!results || tagClosed(PANEL)); }
  // The chart plugin owns mark links and their click handlers. Keep those nodes, but give the child's
  // closed drawing only the producer's fixed geometry. The complete SVG lives off the document until
  // PANEL opens. Older chart files without fixed geometry stay closed in full.
  function veilCharts() {
    Array.prototype.forEach.call(document.querySelectorAll('figure[data-chart]'), function (fig) {
      var sp = fig.pgpSpec, shut = chartPanelClosed();
      if (!sp || fig.pgpPanelState === shut) { return; }
      fig.pgpPanelState = shut;
      if (!fig.pgpFullSVG) { fig.pgpFullSVG = fig.querySelector('svg'); }
      var current = fig.querySelector('svg,.chart-panel-closed');
      if (!shut) {
        if (current && current !== fig.pgpFullSVG) { current.remove(); }
        if (fig.pgpFullSVG && fig.pgpFullSVG.parentNode !== fig) { fig.appendChild(fig.pgpFullSVG); }
        return;
      }
      var replacement = el('p', { 'class': 'dim chart-panel-closed' },
        'This chart stays closed until every result to act on is opened.');
      if (sp.svg_closed && fig.pgpFullSVG) {
        var doc = new DOMParser().parseFromString(sp.svg_closed, 'image/svg+xml');
        if (doc.documentElement.nodeName === 'svg') {
          replacement = document.importNode(doc.documentElement, true);
          Array.prototype.forEach.call(replacement.querySelectorAll('a[data-mark]'), function (a) {
            var id = a.getAttribute('data-mark');
            var live = Array.prototype.find.call(fig.pgpFullSVG.querySelectorAll('a[data-mark]'), function (n) {
              return n.getAttribute('data-mark') === id;
            });
            // Use a clone in the fixed geometry and delegate its click to the plugin's original node.
            if (!live) { a.remove(); }
            else { a.addEventListener('click', function (e) { e.preventDefault(); live.dispatchEvent(new MouseEvent('click', { bubbles: true, cancelable: true })); }); }
          });
        }
      }
      if (current) { current.replaceWith(replacement); } else { fig.appendChild(replacement); }
    });
  }
  // The plugin places figures after person.js draws the slots. Mutation callbacks run before paint,
  // including on the first async charts response; the card event also handles PANEL changing without
  // any mark changing. No chart filename or data value is used to decide whether PANEL is closed.
  new MutationObserver(veilCharts).observe(document.body, { childList: true, subtree: true });
  document.addEventListener('pgp:cards', function () { queueMicrotask(veilCharts); });
  // A card's table without its closed rows: [table or null, how many rows it left out].
  function openTable(c) {
    var tags = veilTags('rows', c.id), t = c.table;
    if (!t || !t.rows) { return [t, 0]; }
    var rows = t.rows.filter(function (r, i) { return !placeClosed(tags[i]); });
    var out = {};
    Object.keys(t).forEach(function (k) { out[k] = t[k]; });
    out.rows = rows;
    return [rows.length ? out : null, t.rows.length - rows.length];
  }
  // A card's facts and drawing without their closed items (a fact, a CHIP call or a copy-number dot that names a
  // closed result's gene): [facts, widget, how many items they left out].
  function openFacts(c) {
    var ft = veilTags('facts', c.id), facts = (c.facts || []).filter(function (f, i) { return !placeClosed(ft[i]); });
    var w = c.widget, key = w && (w.kind === 'vaf' ? 'calls' : w.kind === 'copies' ? 'genes' : null);
    var wt = veilTags('widget', c.id), shut = (c.facts || []).length - facts.length;
    if (key && wt.length && w[key]) {
      var items = w[key].filter(function (x, i) { return !(wt[i] && tagClosed(wt[i])); });
      shut += w[key].length - items.length;
      var copy = {};
      Object.keys(w).forEach(function (k) { copy[k] = w[k]; });
      copy[key] = items;
      w = copy;
    }
    if (w && placeClosed(((results.veil || {}).widget_panel || {})[c.id])) { w = null; shut += 1; }
    return [facts, w, shut];
  }
  function veiledOf(c) { return (c.findings || []).map(findingById).filter(function (f) { return f && veiled(f); }); }
  // The results to act on this page keeps closed, drawn yet or not.
  function closedActs() { return ((results && results.findings) || []).filter(function (x) { return x.level === 'act' && veiled(x); }); }
  // A card of a fixed panel (results.veil.panels, bin/live_results.py FIXED_PANELS: the same genes or loci for
  // everyone, known in advance): a row or dot of it closed with a result would tell the result's gene by what is
  // missing. So it closes whole while this page keeps any result to act on closed, behind one box with the same
  // words whatever the result, and comes back once the last one is opened.
  function panelShut(c) {
    var p = (results && results.veil && results.veil.panels) || [];
    return closedPage() && p.indexOf(c.id) >= 0 && closedActs().length > 0;
  }
  var PANEL_WORDS = 'This section lists the same genes for everyone, so it stays closed while a result to act on is closed ' +
    'on this page: what it left out would tell which gene. It opens once every such result is opened.';
  var PANEL_LINE = 'The same genes for everyone, closed while a result about a child is closed.';
  // The next closed result in the page's order (the map's Open opens the same one).
  function nextClosed() { var ids = mapClosed(); return ids.length ? findingById(ids[0]) : closedActs()[0]; }
  // A closed finding's handle in the page: its place in the page's list of closed findings, never its id, which
  // names its gene and variant. The id comes back once the finding is opened on purpose (#f/<id> too).
  function closedOrder() {
    var out = [];
    sectionsList().forEach(function (sec) {
      (sec.findings || []).forEach(function (id) { var f = findingById(id); if (f && veiled(f) && out.indexOf(id) < 0) { out.push(id); } });
    });
    ((results && results.findings) || []).forEach(function (f) { if (veiled(f) && out.indexOf(f.id) < 0) { out.push(f.id); } });
    return out;
  }
  function handle(f) { return veiled(f) ? 'closed-' + (closedOrder().indexOf(f.id) + 1) : f.id; }
  // The search box (plugins/search.js) asks the same question for each hit, by the summary key and the anchor
  // of the section it comes from: a hit from a card this page keeps closed is left out. A closed page with
  // no results.json yet cannot tell which cards those are, so it leaves every hit out until it can.
  // Asked with no key, it answers whether this page keeps anything closed at all.
  // A hit tagged with a result to act on (veil, bin/results_index.py) is left out while that result is closed.
  window.PGP.searchHidden = function (key, anchor, tag) {
    if (!closedPage()) { return false; }
    if (!results || key === null) { return true; }
    if (tag && tagClosed(tag)) { return true; }
    return cards().some(function (c) {
      return (c.id === anchor || (c.keys || []).indexOf(key) >= 0) && (veiledOf(c).length > 0 || panelShut(c));
    });
  };

  function card(c) {
    var g = gate(c);
    var closedFs = g === 'show' ? veiledOf(c) : [];
    var art = el('article', { 'class': 'card result g-' + g + (attentionNow(c) ? ' attn' : ''), id: 'card-' + c.id,
      'data-section': c.id, 'data-step': (c.steps || []).join(' '), 'data-state': c.state, 'data-level': c.level || '',
      tabindex: '-1', 'aria-labelledby': 'h-' + c.id });
    var head = el('header', { 'class': 'card-head' });
    var h = el('h3', { id: 'h-' + c.id });
    add(h, el('a', { href: '#s/' + c.id, 'class': 'self' }, c.title));
    add(head, h);
    var chips = el('div', { 'class': 'chips' });
    if (g === 'show') {
      // The card's level, on the one scale: a card the report flags is "worth a look" at least.
      var lv = c.level || (attentionNow(c) ? 'review' : 'info');
      add(chips, U.level(lv, attentionNow(c) && lv === 'info' ? 'Worth a look' : undefined));
    }
    recheckWords(c).forEach(function (t) { add(chips, chip(t, 'chip-recheck', '↻')); });
    add(head, chips);
    add(art, head);
    // The text at reading width; the caveat, the numbers and the links in the margin beside it (one
    // column below 1100 px).
    var body = el('div', { 'class': 'card-body' });
    var text = el('div', { 'class': 'card-text' });
    var margin = el('aside', { 'class': 'card-margin' });
    add(body, text, margin);
    add(art, body);
    if (pending(g) && state && state.run && state.run.state === 'stopped') {
      add(text, el('p', { 'class': 'notassessed' }, 'Not written: the run stopped before ' + stepWords(c) + ' ended.'));
    } else if (pending(g)) {
      add(text, el('p', { 'class': 'state-text' }, (g === 'running' ? 'Being written: ' : 'Waiting for ') + stepWords(c) +
        (g === 'running' ? ' is running now.' : '.')));
    } else if (g === 'notrun') {
      add(text, el('p', { 'class': 'notassessed' }, 'Not assessed: no step of this run made this section (step ' + (c.steps || []).join(', ') + ').'));
    } else if (g === 'failed') {
      add(text, add(el('p', { 'class': 'notassessed' }), el('b', {}, 'Not assessed. '),
        'Step ' + (c.steps || []).join(', ') + ' failed: the step timeline says why, and a rerun of the step fills this section.'));
    } else if (g === 'empty') {
      add(text, el('p', { 'class': 'notassessed' }, 'Not assessed: its step ended without a result for this section.'));
    } else if (panelShut(c)) {
      // A fixed-panel card says the same words whichever result is closed, its own findings' too: so it comes first.
      // One button opens one result, the next closed one in the page's order, as the map's Open does.
      add(text, U.veil(PANEL_WORDS, null, function () { var f = nextClosed(); if (f) { openFinding(f); } redraw(true); }));
    } else if (closedFs.length) {
      // One button opens one result: a card with two results to act on asks twice.
      add(text, U.veil((closedFs.length === 1 ? 'This section holds a result' : 'This section holds ' + closedFs.length + ' results') +
        ' to act on in adult life. Its rows stay closed until ' + (closedFs.length === 1 ? 'it is' : 'each is') + ' opened.', null,
        function () { openFinding(closedFs[0]); redraw(true); var x = document.getElementById('f-' + closedFs[0].id) || document.getElementById('card-' + c.id); if (x) { x.scrollIntoView({ block: 'start' }); } }));
    } else {
      (c.notes || []).forEach(function (n, i) { if (!placeClosed(veilTags('notes', c.id)[i])) { add(text, el('p', { 'class': 'note-stale' }, n)); } });
      if (c.headline && !headlineShut(c)) { add(text, el('p', { 'class': 'headline' }, c.headline)); }
      var of = openFacts(c);
      add(margin, U.facts(of[0]));
      add(text, widget(of[1]));
      Object.keys((results && results.charts) || {}).forEach(function (chart) {
        if (results.charts[chart] === c.id) { add(text, el('div', { 'class': 'chart-slot', 'data-chart-slot': chart })); }
      });
      if (c.image) { add(text, el('img', { 'class': 'plot', src: c.image.src, alt: c.image.alt })); }
      var ot = openTable(c), tb = U.table(ot[0]);
      var lt = veilTags('lines', c.id);
      var lines = (c.lines || []).filter(function (l, i) { return !placeClosed(lt[i]); });
      var shut = ot[1] + (c.lines || []).length - lines.length + (c.headline && headlineShut(c) ? 1 : 0), items = of[2];
      if (tb) { add(text, tb); } else if (c.empty_text && !ot[1]) { add(text, el('p', { 'class': 'empty' }, c.empty_text)); }
      if (shut || items) {
        // The same words however many rows or items close, and whichever result they belong to: a count would
        // tell one closed result from another.
        add(text, el('p', { 'class': 'dim veil-rows' }, 'Some of its rows or numbers stay closed with a result to act on, until that result is opened.'));
      }
      add(text, lineList(lines));
      if (c.hidden) {
        // APOE and the like: behind their warning, opened on purpose.
        // Its fold is the map's APOE switch too (map.js): opening it here opens APOE on the map, and back.
        var d = el('details', { 'class': 'behind', 'data-gate': 'apoe' });
        add(d, el('summary', {}, c.hidden.summary), el('p', {}, c.hidden.text));
        d.addEventListener('toggle', function () { if (d.open !== !!shown.apoe) { setShown('apoe', d.open); } });
        add(text, d);
      }
      if (c.explore && c.explore.length) {
        var p = el('p', { 'class': 'explore' }, 'Explore the branch (opens a website; the address carries the haplogroup): ');
        c.explore.forEach(function (x, i) {
          if (i) { add(p, ', '); }
          add(p, el('a', { href: x.url, rel: 'noreferrer noopener', target: '_blank' }, x.name));
        });
        add(text, p);
      }
    }
    add(margin, U.caveat(c.caveat));
    var foot = el('p', { 'class': 'card-links' });
    add(foot, el('a', { href: 'report.html#' + encodeURIComponent(c.anchor), 'class': 'report' }, 'Open in the report'));
    var vtags = veilTags('viewer', c.id);
    var reads = (c.viewer || []).filter(function (x, i) { return !placeClosed(vtags[i]); });
    if (g === 'show' && reads.length && !panelShut(c)) {
      add(foot, el('a', { href: 'viewer/#' + reads[0], 'class': 'viewer' }, 'See the reads at ' + reads[0]));
    }
    (g === 'show' ? c.links || [] : []).forEach(function (f) {
      add(foot, el('a', { href: f.split('/').map(encodeURIComponent).join('/') }, 'Open ' + f));
    });
    add(margin, foot);
    return art;
  }
  // A card's lines. Up to LINES_SHOWN show as they are. A longer list keeps the lines that say something
  // in view and folds the rest (an allele absent, a call not determined) behind one line that counts them;
  // every line stays in the page, so a search or a reader still finds it.
  var LINES_SHOWN = 6;
  var QUIET_LINE = /\b(absent|not determined|not typed|not called|not found)\b/i;
  function lineList(lines) {
    if (!lines || !lines.length) { return null; }
    var keep = lines, rest = [];
    if (lines.length > LINES_SHOWN) {
      keep = lines.filter(function (l) { return !QUIET_LINE.test(l); });
      rest = lines.filter(function (l) { return QUIET_LINE.test(l); });
      if (keep.length > LINES_SHOWN) { rest = keep.slice(LINES_SHOWN).concat(rest); keep = keep.slice(0, LINES_SHOWN); }
    }
    var box = el('div', { 'class': 'lines-box' });
    var ul = el('ul', { 'class': 'lines' });
    keep.forEach(function (l) { add(ul, el('li', {}, l)); });
    if (keep.length) { add(box, ul); }
    if (rest.length) {
      var quiet = rest.filter(function (l) { return QUIET_LINE.test(l); }).length;
      var d = el('details', { 'class': 'lines-more' });
      add(d, el('summary', {}, (keep.length ? 'The other ' + rest.length + ' of ' + lines.length + ' lines' : 'All ' + lines.length + ' lines') +
        (quiet === rest.length ? ': absent or not determined' : quiet ? ' (' + quiet + ' absent or not determined)' : '')));
      var ul2 = el('ul', { 'class': 'lines' });
      rest.forEach(function (l) { add(ul2, el('li', {}, l)); });
      add(d, ul2);
      add(box, d);
    }
    return box;
  }

  // --- a finding in its section ------------------------------------------------------------------------
  function findingHref(f) { return 'finding/' + encodeURIComponent(f.id) + '/'; }
  function newAt(f) {
    if (!running()) { return null; }
    var by = stepsBy(), t = null;
    (f.steps || []).forEach(function (s) { if (by[s] && by[s].ended_utc) { t = by[s].ended_utc; } });
    return t;
  }
  function note(head, body, cls) {
    var n = el('div', { 'class': 'note' + (cls ? ' ' + cls : '') });
    add(n, el('span', { 'class': 'nh' }, head));
    add(n, body);
    return n;
  }
  function evidence(rows) {
    if (!rows || !rows.length) { return null; }
    var dl = el('dl', {});
    rows.forEach(function (r) { add(dl, el('dt', {}, r[0]), el('dd', {}, r[1])); });
    return dl;
  }
  // key: the finding's id, so a redraw keeps this fold open on its own finding (foldKeys).
  function caveatNote(list, key) {
    if (!list || !list.length) { return null; }
    var box = el('div', {});
    var lead = U.firstSentence(list[0]);
    add(box, el('span', {}, lead + ' '));
    var d = el('details', key ? { 'data-card': 'cav-' + key } : {});
    add(d, el('summary', { 'class': 'toggle' }, 'More'));
    list.forEach(function (t) { add(d, el('p', {}, t)); });
    add(box, d);
    return note('Caveat', box, 'caveat-note');
  }
  // A finding that is not itself a result to act on but names one (a symptom ranking can rank the gene: results.veil.findings)
  // stays closed while that result is, and opening it opens that result.
  function findingTag(f) {
    if (f.level === 'act') { return null; }  // its own switch; never an unclassified secondary place
    var t = ((results && results.veil && results.veil.findings) || {})[f.id];
    if (t === undefined) { t = ((results && results.findings) || []).some(function (x) { return x.level === 'act'; }) ? [PANEL] : null; }
    t = (Array.isArray(t) ? t : t ? [t] : []).filter(function (x) { return x !== f.id; });
    return t.length ? t : null;
  }
  function veiled(f) {
    if (!closedPage()) { return false; }
    if (f.level === 'act') { return !shown[f.id]; }
    var t = findingTag(f);
    return !!(t && tagClosed(t));
  }
  // Opening a finding opens every result its tag names. One tagged PANEL (a fixed-panel card's own) stays closed while
  // any result to act on is: opening it opens the next closed result in the page's order too, as its card's box does.
  function openFinding(f) {
    var t = findingTag(f) || [];
    keepHash(function () {
      shown[f.id] = true;
      t.forEach(function (x) { if (x !== PANEL) { shown[x] = true; } });
      var n = t.indexOf(PANEL) >= 0 ? nextClosed() : null;
      if (n) { shown[n.id] = true; }
    });
  }
  // The address never holds the id of a finding this page can close, closed or opened: #f/closed-<n> while it is
  // closed (its handle) and #f/opened-<n> once opened on purpose (its place among the opened ones), so no history
  // entry keeps the id. Any other #f/<x> on a page that keeps a result closed, a closed result's id or no id at all,
  // lands the same way, on #f/closed-1: the page never picks a closed card by its id, so a guessed id tells nothing.
  // Opening or closing a result moves the handles of the others: the address follows the finding it named. The
  // finding behind a handle in the address is kept here only (hashId); a new address from the reader forgets it.
  var hashId = null;
  window.addEventListener('hashchange', function () { hashId = null; });
  function closable(f) { return closedPage() && (f.level === 'act' || !!findingTag(f)); }
  function openedOrder() {
    var out = [];
    sectionsList().forEach(function (sec) {
      (sec.findings || []).forEach(function (id) { var f = findingById(id); if (f && closable(f) && !veiled(f) && out.indexOf(id) < 0) { out.push(id); } });
    });
    ((results && results.findings) || []).forEach(function (f) { if (closable(f) && !veiled(f) && out.indexOf(f.id) < 0) { out.push(f.id); } });
    return out;
  }
  function addr(f) { return !closable(f) ? f.id : veiled(f) ? handle(f) : 'opened-' + (openedOrder().indexOf(f.id) + 1); }
  function hashFinding() {
    var h;
    try { h = decodeURIComponent((window.location.hash || '').slice(1)); } catch (e) { return null; }
    if (h.indexOf('f/') !== 0) { return null; }
    var x = h.slice(2), m = /^(closed|opened)-(\d+)$/.exec(x);
    var f = (m && hashId && findingById(hashId)) ||
      (m ? findingById((m[1] === 'closed' ? closedOrder() : openedOrder())[Number(m[2]) - 1]) : findingById(x));
    if (f && !(veiled(f) && !m)) { return f; }
    return closedOrder().length ? findingById(closedOrder()[0]) : null;
  }
  function syncHash() {
    var f = hashFinding();
    if (!f || !window.history || !window.history.replaceState) { return; }
    hashId = f.id;
    var to = '#f/' + addr(f);
    if (window.location.hash !== to) { window.history.replaceState(null, '', to); }
  }
  function keepHash(change) { var at = hashFinding(); if (at) { hashId = at.id; } change(); syncHash(); }
  // The genome map's veil (map.js): 'closed' while this page keeps the result closed, 'opened' once it was
  // opened on purpose, 'open' when this page never closes it. A gated result (APOE, behind its warning) is
  // closed on every page, an adult's too, until opened; a tag naming no finding of these results stays closed.
  function mapVeil(id, gated) {
    if (gated) { return shown[id] ? 'opened' : 'closed'; }
    if (!closedPage()) { return 'open'; }
    var f = findingById(id);
    if (!f || veiled(f)) { return 'closed'; }
    return f.level === 'act' || findingTag(f) ? 'opened' : 'open';
  }
  // One switch for a result, from the map's chip or the card's "Show the result" (openFinding): opening it in
  // one place opens it in both. Closing it closes the findings that were opened with it.
  function setShown(id, on) {
    var f = findingById(id);
    if (on) { if (f) { openFinding(f); } else { shown[id] = true; } } else {
      keepHash(function () {
        delete shown[id];
        ((results && results.findings) || []).forEach(function (x) { if ((findingTag(x) || []).indexOf(id) >= 0) { delete shown[x.id]; } });
      });
    }
    redraw(true);
  }
  // Whether the page shows a finding of the map at this moment: its card's steps are done (a replay draws
  // the run's end results from its first frame). A row with no finding of these results (the map's own APOE
  // row) follows the card that holds its fold (the traits card): a replay before that card shows has no APOE.
  // A finding this page keeps closed is not shown either, whatever row tag it shares with another result.
  function mapNow(id) {
    var f = findingById(id);
    if (!f) { var t = cards().filter(function (c) { return c.hidden; })[0]; return !!(t && gate(t) === 'show'); }
    var c = cardById(f.card);
    return !!(c && gate(c) === 'show') && !veiled(f);
  }
  // The results to act on this page keeps closed at this moment, in the page's order: the map's chip counts
  // them, whether or not they have a place on the map.
  function mapClosed() {
    return findingsNow().filter(function (f) { return f.level === 'act' && veiled(f); }).map(function (f) { return f.id; });
  }
  // The title of a result once the page shows it, for the map's "Opened in its section" chip; never a closed one's.
  // The address of a result the map opened, for its "Opened in its section" link: never its id (addr()).
  function mapHref(id) { var f = findingById(id); return f ? addr(f) : id; }
  function mapTitle(id) { var f = findingById(id); return f && !veiled(f) ? f.title : ''; }
  // Whether a row of the map is a finding of these results: then the page's own veil alone decides it (map.js).
  function mapKnown(id) { return !!findingById(id); }
  // Whether a step of the run has ended at this moment (a replay before it has not): the map's layers follow it.
  function stepDone(step) { var s = stepsBy()[step]; return !s || s.state === 'done'; }
  // Without the map (a browser with no DecompressionStream, or a map the server could not build): the findings
  // with a place, as a table, each through this page's veil first.
  function mapFallback() {
    return findingsNow().filter(function (f) { return f.viewer && !veiled(f); }).map(function (f) {
      return [U.LEVEL_WORDS[f.level] || f.level, f.title, f.viewer];
    });
  }
  function findingArticle(f, big) {
    var row = el('div', { 'class': 'row' });
    var art = el('div', { 'class': 'finding', id: 'f-' + handle(f), 'data-finding': handle(f), 'data-level': f.level });
    var at = newAt(f);
    if (at) { art.classList.add('new'); add(art, el('span', { 'class': 'newband' }, 'New at ' + U.hm(at))); add(art, el('br', {})); }
    add(art, U.findingLevel(f.level, f.level_words));
    if (veiled(f)) {
      add(art, U.veil('A result to act on in adult life is usually not looked at before the child can choose. It stays closed until you open it.',
        null, function () { openFinding(f); redraw(true); var x = document.getElementById('f-' + f.id); if (x) { x.scrollIntoView({ block: 'start' }); } }));
      add(row, art, el('aside', {}));
      return row;
    }
    var h = el('h3', { 'class': big ? 'h1' : 'h3' });
    add(h, el('a', { href: findingHref(f) }, f.title));
    add(art, h);
    if (f.detail) { add(art, el('div', { 'class': 'gene' }, f.detail)); }
    add(art, el('p', { 'class': 'body-serif' }, f.text));
    if (big && f.plan && f.plan.length) {
      add(art, el('div', { 'class': 'todo-h' }, 'What to do'));
      var ol = el('ol', { 'class': 'todo' });
      f.plan.forEach(function (t) { add(ol, el('li', {}, t)); });
      add(art, ol);
    }
    add(art, el('a', { 'class': 'readmore', href: findingHref(f) }, 'The full evidence and what to do →'));
    var side = el('aside', {});
    var ev = evidence(f.evidence);
    if (ev) { add(side, note('Evidence', ev)); }
    add(side, caveatNote(f.caveat, f.id));
    var rel = (f.related || []).map(findingById).filter(function (r) { return r && !veiled(r); });
    if (rel.length) {
      var box = el('div', {});
      rel.slice(0, 4).forEach(function (r, i) { if (i) { add(box, el('br', {})); } add(box, el('a', { href: '#f/' + addr(r) }, r.title)); });
      add(side, note('Linked', box, 'linked-note'));
    }
    add(row, art, side);
    return row;
  }
  function findingLine(f) {
    var li = el('li', { id: 'f-' + handle(f), 'data-finding': handle(f), 'data-level': f.level });
    if (veiled(f)) {
      // A finding that names a closed result (a symptom ranking of its gene): closed with it, opened with it.
      var b = el('button', { type: 'button', 'class': 'btn' }, 'Show the result');
      b.addEventListener('click', function () { openFinding(f); redraw(true); var x = document.getElementById('f-' + f.id); if (x) { x.scrollIntoView({ block: 'start' }); } });
      add(li, el('span', { 'class': 't' }, 'A result about a child, closed until opened.'), U.findingLevel(f.level, f.level_words), b);
      return li;
    }
    add(li, el('a', { 'class': 't', href: findingHref(f) }, f.title), U.findingLevel(f.level, f.level_words));
    if (f.detail) { add(li, el('span', { 'class': 'gene' }, f.detail)); }
    add(li, el('span', { 'class': 'x' }, f.text));
    return li;
  }

  // --- the sections ---------------------------------------------------------------------------------
  // A section's stamp: written (with the time its last step ended), being written (a step runs), waiting.
  // The chapters whose steps can still add a finding to sections 1 to 4 while they run.
  var FEEDS = { clinical: true, carrier: true, medicines: true };
  function feeding(sec) {
    return sec.n <= 4 ? cards().filter(function (c) { return FEEDS[c.chapter] && sec.cards.indexOf(c.id) < 0; }) : [];
  }
  function stamp(sec) {
    var cs = sec.cards.map(cardById).filter(Boolean).concat(feeding(sec));
    if (!cs.length) { return { k: 'written', t: 'Written' }; }
    var gs = cs.map(gate);
    if (gs.some(function (g) { return g === 'running'; })) { return { k: 'writing', t: 'Being written' }; }
    if (gs.some(function (g) { return g === 'waiting'; })) {
      return gs.some(function (g) { return g === 'show'; }) ? { k: 'writing', t: 'Half written' } : { k: 'waiting', t: 'Waiting' };
    }
    var by = stepsBy(), t = null;
    cs.forEach(function (c) { (c.steps || []).forEach(function (s) { if (by[s] && by[s].ended_utc && (!t || by[s].ended_utc > t)) { t = by[s].ended_utc; } }); });
    return { k: 'written', t: 'Written' + (t && running() ? ' ' + U.hm(t) : '') };
  }
  function sectionHead(sec) {
    var head = el('div', { 'class': 'sechead' });
    add(head, el('span', { 'class': 'secnum' }, String(sec.n)), el('h2', { 'class': 'h2' }, sec.title));
    if (running()) { var st = stamp(sec); add(head, el('span', { 'class': 'stamp ' + st.k }, st.t)); }
    return head;
  }
  function waitBlock(sec) {
    var cs = sec.cards.map(cardById).filter(function (c) { return c && pending(gate(c)); });
    if (!cs.length) { return null; }
    var box = el('div', {});
    add(box, add(el('div', { 'class': 'proof', 'aria-hidden': 'true' }), el('i', {}), el('i', {}), el('i', {})));
    var names = cs.map(function (c) { return c.title + (gate(c) === 'running' ? ' (running)' : ''); });
    add(box, el('p', { 'class': 'waitnote' }, 'This section stays open until its steps finish: ' + names.join(', ') +
      '. Until then this page says nothing either way.'));
    return box;
  }
  // What the clinical steps of section 8 looked at, named in an empty section 1 to 4.
  function assessedNames() {
    var s8 = sectionsList().filter(function (s) { return s.n === 8; })[0];
    return s8 ? s8.cards.map(cardById).filter(function (c) { return c && gate(c) === 'show'; }).map(function (c) { return c.title; }) : [];
  }
  function sectionBody(sec) {
    var body = el('div', {});
    var fs = sec.findings.map(findingById).filter(function (f) { var c = f && cardById(f.card); return c && gate(c) === 'show'; });
    var cs = sec.cards.map(cardById).filter(Boolean);
    var big = sec.n === 1 || (sec.n === 4 && fs.length <= 2);
    if (sec.n <= 4 || sec.n === 5) {
      if (fs.length) {
        if (big || sec.n === 1 || fs.length <= 3) { fs.forEach(function (f) { add(body, findingArticle(f, big)); }); }
        else { var ul = el('ul', { 'class': 'list-findings' }); fs.forEach(function (f) { add(ul, findingLine(f)); }); add(body, ul); }
      } else if (sec.n !== 5) {
        var w = waitBlock(sec);
        var pend = feeding(sec).some(function (c) { return pending(gate(c)); });
        if (w) { add(body, w); } else if (pend) {
          add(body, el('p', { 'class': 'waitnote' }, 'This section stays open until the steps that can add a finding here finish. Until then this page says nothing either way.'));
        } else {
          var e = el('p', { 'class': 'empty-sec' }, NOT_ALL_CLEAR);
          add(body, e);
          var names = sec.n === 1 || sec.n === 4 ? assessedNames() : cs.filter(function (c) { return gate(c) === 'show'; }).map(function (c) { return c.title; });
          if (names.length) { add(body, el('p', { 'class': 'soft' }, 'Assessed: ' + names.join(', ') + '. Section 8 lists what each looked at, section 9 what no step sees.')); }
        }
      }
    }
    if (sec.n === 8) {
      var nf = el('ul', { 'class': 'nf' });
      cs.forEach(function (c) {
        var li = el('li', {});
        var d = el('details', { 'data-card': c.id });
        var sum = el('summary', {});
        var w2 = el('span', { 'class': 'what' }, cardLine(c));
        if (c.caveat && c.caveat.length) { add(w2, el('span', { 'class': 'cav' }, 'Does not see: ' + U.firstSentence(c.caveat[0]))); }
        add(sum, el('span', { 'class': 'nfname' }, c.title), w2);
        add(d, sum, card(c));
        add(li, d);
        add(nf, li);
      });
      add(body, nf);
      return body;
    }
    if ((sec.n === 9 || sec.n === 10) && tagClosed(PANEL)) {
      add(body, U.veil('This part can contain data text, so it stays closed until every result to act on is opened.',
        null, function () { var f = nextClosed(); if (f) { openFinding(f); } redraw(true); }));
      return body;
    }
    if (sec.n === 9) {
      var lists = [['No step looks at these', (results && results.not_assessed) || []], ['Steps that did not run', (results && results.not_run) || []]];
      lists.forEach(function (x) {
        add(body, el('h3', { 'class': 'h3' }, x[0]));
        var ul2 = el('ul', { 'class': 'cols2' });
        (x[1].length ? x[1] : ['None listed.']).forEach(function (t) { add(ul2, el('li', {}, t)); });
        add(body, ul2);
      });
      return body;
    }
    if (sec.n === 10) {
      var dbs = (results && results.databases) || {};
      var keys = Object.keys(dbs);
      if (keys.length) {
        var t = el('table', { 'class': 'sources' });
        add(t, add(el('tr', {}), el('th', {}, 'Source'), el('th', {}, 'Version or date')));
        keys.forEach(function (k) { add(t, add(el('tr', {}), el('td', {}, k.replace(/_/g, ' ')), el('td', {}, dbs[k]))); });
        add(body, t);
      }
      var rc = state && state.plugins && state.plugins.recheck;
      if (rc && rc.banner) { add(body, el('p', { 'class': 'note-stale' }, rc.banner)); }
      var fb = el('div', { 'class': 'filebtns' });
      add(fb, el('a', { 'class': 'btn', href: 'report.html' }, 'The full report (HTML)'));
      add(fb, el('a', { 'class': 'btn ghost', href: 'viewer/' }, 'Genome viewer'));
      add(body, fb);
    }
    cs.forEach(function (c) { add(body, phoneFold(c, sec) ? cardFold(c) : card(c)); });
    return body;
  }

  // --- the phone: what matters open, the rest one tap away -------------------------------------------
  // Under 700 px the findings of sections 1 to 4 stay open; the cards behind them, every information
  // card, and the bodies of sections 5 to 10 start closed, each behind one line that says what it holds.
  // A link to a card, a section or a finding opens the folds around it (route). Wider screens draw
  // every card open, as before.
  var phoneQuery = window.matchMedia ? window.matchMedia('(max-width: 700px)') : null;
  function phone() { return !!(phoneQuery && phoneQuery.matches); }
  var LOUD = { act: true, medicine: true, carrier: true, review: true, trust: true };
  // One line for a card: what it found, or why it has nothing to show yet.
  function cardLine(c) {
    var g = gate(c);
    return g === 'show' && veiledOf(c).length ? 'A result about a child, closed until opened.'
      : g === 'show' && panelShut(c) ? PANEL_LINE
      : g === 'show' ? ((!headlineShut(c) && c.headline) || (openFacts(c)[0][0] ? openFacts(c)[0][0][0] + ': ' + openFacts(c)[0][0][1] : NOT_ALL_CLEAR))
      : pending(g) ? 'Waiting for ' + stepWords(c) + '.' : g === 'failed' ? 'Not assessed: its step failed.' : 'Not assessed in this run.';
  }
  // A card the report flags that has no finding of its own: on a phone it is the only place that flag
  // shows, so it stays open, in section 5 too.
  function flaggedAlone(c) { return attentionNow(c) && !(c.findings || []).some(function (f) { return findingById(f); }); }
  // A card folds on a phone when it only informs, or when its findings already stand open in sections 1
  // to 4 (the card is their evidence). A card the report flags with no finding of its own stays open.
  function phoneFold(c, sec) {
    if (!phone() || sec.n === 8) { return false; }
    if (sec.n >= 5) { return !(sec.n === 5 && flaggedAlone(c)); }
    if (!(LOUD[c.level] || attentionNow(c))) { return true; }
    var open = {};
    sectionsList().forEach(function (s) { if (s.n <= 4) { s.findings.forEach(function (f) { open[f] = true; }); } });
    var fs = (c.findings || []).filter(function (f) { return findingById(f); });
    return fs.length > 0 && fs.every(function (f) { return open[f]; });
  }
  function foldSummary(name, line) {
    var sum = el('summary', {});
    if (name) { add(sum, el('span', { 'class': 'fname' }, name)); }
    add(sum, el('span', { 'class': 'what' }, line));
    return sum;
  }
  function cardFold(c) {
    var d = el('details', { 'class': 'pfold', 'data-card': 'fold-' + c.id });
    add(d, foldSummary(c.title, cardLine(c)), card(c));
    return d;
  }
  // Sections 5 to 10 fold whole on a phone, unless section 5 holds a finding (a data-quality problem) or
  // a card the report flags with no finding of its own.
  function sectionFolds(sec) {
    if (!phone() || sec.n < 5) { return false; }
    return !(sec.n === 5 && (sec.findings.some(function (id) { var f = findingById(id), c = f && cardById(f.card); return c && gate(c) === 'show'; }) ||
      sec.cards.map(cardById).some(function (c) { return c && flaggedAlone(c); })));
  }
  function sectionLine(sec) {
    if (sec.n === 9 && tagClosed(PANEL)) { return 'Closed until every result to act on is opened.'; }
    if (sec.n === 9) {
      var na = ((results && results.not_assessed) || []).length, nr = ((results && results.not_run) || []).length;
      return na + (na === 1 ? ' thing' : ' things') + ' no step looks at, ' + nr + (nr === 1 ? ' step' : ' steps') + ' that did not run';
    }
    var names = sec.cards.map(cardById).filter(Boolean).map(function (c) { return c.title; });
    if (sec.n === 10) { names = ['Database versions', 'the full report'].concat(names); }
    return names.length ? names.join(', ') : NOT_ALL_CLEAR;
  }
  function drawSections() {
    panel.textContent = '';
    if (results && results.note && !(results.chapters || []).length) {
      add(panel, el('p', { 'class': 'empty' }, results.note));
      return;
    }
    sectionsList().forEach(function (sec) {
      var s = el('section', { 'class': 'sec', id: 's' + sec.n, 'data-journal': sec.key, 'aria-labelledby': 'sh' + sec.n });
      var head = sectionHead(sec);
      head.querySelector('h2').id = 'sh' + sec.n;
      add(s, head, el('div', { 'class': 'secq' }, sec.question));
      if (sectionFolds(sec)) {
        var d = el('details', { 'class': 'pfold secfold', 'data-card': 'sec-' + sec.n });
        add(d, foldSummary('', sectionLine(sec)), sectionBody(sec));
        add(s, d);
      } else { add(s, sectionBody(sec)); }
      add(panel, s);
    });
    root.setAttribute('data-pgp-results', String(panel.querySelectorAll('article.result').length));
    root.setAttribute('data-pgp-findings', String(panel.querySelectorAll('[data-finding]').length));
    document.dispatchEvent(new CustomEvent('pgp:cards'));
  }

  // --- the lead, the counts by level, the chromosome map ---------------------------------------------
  function counts() {
    var n = { act: 0, medicine: 0, review: 0, carrier: 0, trust: 0, info: 0 };
    findingsNow().forEach(function (f) { n[f.level] += 1; });
    n.info = cards().filter(function (c) { return gate(c) === 'show' && c.level === 'info'; }).length;
    return n;
  }
  // A segment of the lead that names a closed result, without the names: the list after a colon, the genes in
  // brackets, or what follows "of".
  function closedLead(t) {
    var w = 'closed on this page until opened on purpose, in section 1';
    if (/^\s*:/.test(t)) { return ': ' + w + '. '; }
    if (/\(.*\)/.test(t)) { return t.replace(/\(.*\)/, '(' + w + ')'); }
    if (/ of /.test(t)) { return t.replace(/ of .*$/, ' result, ' + w); }
    return w + '. ';
  }
  function drawLead() {
    var slot = document.getElementById('lead-slot');
    slot.textContent = '';
    if (!results || !(results.chapters || []).length) { return; }
    var row = el('div', { 'class': 'row' });
    var main = el('div', {});
    var p = el('p', { 'class': 'lede dropcap', id: 'lead' });
    var closed = (results.findings || []).some(veiled);
    var leadTags = (results.veil && results.veil.lead) || [];
    var parts = (results.lead || []).map(function (x, i) {
      if (x[2] === 'act' && closed) { return [': closed on this page until opened on purpose, in section 1. ', false]; }
      return leadTags[i] && tagClosed(leadTags[i]) ? [closedLead(x[0]), x[1]] : x;
    });
    if (replay && running()) {
      // A replay holds the results of the finished run: before its end the lead would name findings the
      // sections do not show yet.
      parts = [['This report is being written. ', true], ['The lead paragraph is written when the run ends; the sections below fill in as their steps finish.', false]];
    } else if (running()) {
      parts = [['This report is being written. ', true], ['What follows is what has been found so far; a section with steps still to run says so. ', false]]
        .concat(parts.filter(function (x) { return !/^Mean depth/.test(x[0]); }));
    }
    parts.forEach(function (x) { add(p, x[1] ? el('b', {}, x[0]) : x[0]); });
    add(main, p);
    var n = counts();
    var fs = findingsNow();
    var led = el('div', { 'class': 'ledger', id: 'ledger' });
    [['act', 1], ['medicine', 2], ['review', 4], ['carrier', 3], ['info', 6]].forEach(function (x) {
      var a = el('a', { href: '#s' + x[1], 'data-level': x[0] });
      add(a, U.level(x[0], { act: 'Take to a doctor', medicine: 'Prescriber', review: 'Worth a look', carrier: 'Family planning', info: 'Information' }[x[0]]));
      add(a, el('span', { 'class': 'n' }, String(n[x[0]])));
      var words = x[0] === 'info' ? 'Origins, traits, scores' : x[0] === 'act' && closed && n.act ? 'Closed until opened' :
        fs.filter(function (f) { return f.level === x[0] && !veiled(f); }).slice(0, 2).map(function (f) { return f.gene || f.title; }).join(', ') ||
        (fs.some(function (f) { return f.level === x[0]; }) ? 'Closed until opened' : 'None in what was assessed');
      add(a, el('span', { 'class': 'w' }, words));
      add(led, a);
    });
    add(main, led);
    var trust = el('div', { 'class': 'trustline' });
    var problems = fs.filter(function (f) { return f.level === 'trust'; });
    var openProblems = problems.filter(function (f) { return !veiled(f); });
    add(trust, U.level('trust', 'Data quality', problems.length ? 'problem' : ''));
    var qc = cardById('qc');
    if (qc && gate(qc) === 'show') {
      openFacts(qc)[0].filter(function (f) { return /depth|contamination|sex from/i.test(f[0]); }).slice(0, 3).forEach(function (f, i) {
        if (i) { add(trust, el('span', { 'class': 'sep' }, '·')); }
        add(trust, el('span', { title: f[0] + ': ' + f[1] }, String(f[0]).replace(/\s*\([^)]*\)/g, '') + ': ' + String(f[1]).split(/[.;(]\s/)[0].replace(/[.;]$/, '').trim()));
      });
    } else { add(trust, el('span', {}, 'Waiting for the quality steps.')); }
    openProblems.forEach(function (f) { add(trust, el('span', { 'class': 'sep' }, '·'), el('a', { href: '#f/' + addr(f) }, f.title)); });
    add(trust, el('span', { 'class': 'sep' }, '·'), el('a', { href: '#s5' }, 'Can I trust this?'));
    add(main, trust);
    var side = el('aside', {});
    add(side, note('How to read this', el('span', {}, 'Sections run from what matters most to what is only interesting. A section with nothing to report still says what it looked at and what it cannot see.'), 'howto'));
    add(side, note('A research pipeline', el('span', {}, 'Not a clinical test: anything to act on is confirmed by a certified laboratory first.')));
    add(row, main, side);
    add(slot, row);
    // Figure 1: the genome map (map.js), one figure kept across redraws so its camera stays where the
    // reader left it. It asks this page's veil for each closed row, and opens or closes a result through it.
    if (window.PGPMap) {
      if (!genomeMap) {
        genomeMap = window.PGPMap.create({ url: 'map.json', veil: mapVeil, now: mapNow, known: mapKnown, closed: mapClosed, title: mapTitle, href: mapHref, stepDone: stepDone, open: function (id) { setShown(id, true); },
          close: function (id) { setShown(id, false); }, fallback: mapFallback });
      }
      add(slot, genomeMap.node);
      genomeMap.refresh(resultsKey);
    }
  }

  // --- the table of contents, the phone tab bar, the cover ------------------------------------------
  function drawToc() {
    rail.textContent = '';
    if (!sectionsList().length) { return; }
    add(rail, el('div', { 'class': 'group' }, running() ? 'The report, section by section' : 'In this report'));
    var ol = el('ol', {});
    var n = counts();
    var lvOf = { 1: 'act', 2: 'medicine', 3: 'carrier', 4: 'review' };
    sectionsList().forEach(function (sec) {
      var a = el('a', { href: '#s' + sec.n, 'data-toc': sec.n });
      add(a, el('span', { 'class': 'tn' }, String(sec.n)), el('span', {}, sec.title));
      var lv = lvOf[sec.n];
      if (lv && n[lv]) { add(a, U.level(lv, String(n[lv]), 'tc')); }
      else if (sec.n === 5) { add(a, U.level('trust', '', 'tc' + (n.trust ? ' problem' : ''))); }
      else if (sec.n === 8) { add(a, el('span', { 'class': 'tc muted' }, String(sec.cards.length))); }
      if (running()) { var st = stamp(sec); add(a, el('span', { 'class': 'tsub ' + st.k }, st.t.toLowerCase())); }
      add(ol, add(el('li', {}), a));
    });
    add(rail, ol);
    add(rail, el('div', { 'class': 'group' }, 'How the run went'));
    add(rail, add(el('div', { 'class': 'people' }), add(el('a', { href: '#progress' }), el('span', {}, 'Step timeline'), el('span', { 'class': 'muted' }, '→'))));
    if (family) {
      add(rail, el('div', { 'class': 'group' }, 'The family'));
      var ppl = el('div', { 'class': 'people' });
      add(ppl, add(el('a', { href: family.url }), el('span', {}, 'Family report'), el('span', { 'class': 'muted' }, '→')));
      (family.data.state.members || []).forEach(function (m) {
        if (m.sample === personName()) { return; }
        var u = memberUrl(m.sample);
        var a = el(u ? 'a' : 'span', u ? { href: u } : { 'class': 'muted' });
        add(a, el('span', { 'class': 'pc' + (m.sequenced ? '' : ' nogenome') }, m.sample));
        var r = family.data.runs && family.data.runs[m.sample];
        if (r && r.attention) { add(a, U.level('review', String(r.attention), 'tc')); }
        add(ppl, a);
      });
      add(rail, ppl);
    }
  }
  function drawTabbar() {
    var bar = document.getElementById('tabbar');
    if (!bar) { return; }
    bar.textContent = '';
    var n = counts();
    [[1, 'To act', 'act'], [2, 'Medicines', 'medicine'], [3, 'Family', 'carrier'], [5, 'Trust', 'trust'], [8, 'More', null]].forEach(function (x) {
      var a = el('a', { href: '#s' + x[0] });
      add(a, el('span', { 'class': 'ti' }, String(x[0])));
      add(a, el('span', {}, x[1] + (x[2] && n[x[2]] ? ' · ' + n[x[2]] : '')));
      add(bar, a);
    });
  }
  function header() {
    var p = (results && results.person) || {};
    var planned = ((state && state.steps) || []).filter(function (s) { return s.state !== 'not planned'; });
    var done = planned.filter(function (s) { return s.state === 'done'; }).length;
    var slot = document.getElementById('avatar-slot');
    slot.textContent = '';
    add(slot, U.avatar(p.inferred_sex || (state && state.run && state.run.declared_sex), true, planned.length ? done / planned.length : 0));
    var run = (state && state.run) || {};
    var kicker = document.getElementById('kicker');
    kicker.textContent = (running() ? 'Being written · ' : 'Personal genome report') + (family ? (running() ? '' : ' · ') + 'Family ' + family.name : '');
    var line = document.getElementById('person-line');
    line.textContent = '';
    var bits = [];
    if (p.inferred_sex || p.declared_sex) {
      bits.push(el('span', {}, (p.inferred_sex ? p.inferred_sex + ' from the reads' : 'sex not measured') +
        (p.declared_sex ? (p.declared_sex === p.inferred_sex ? ', as declared' : ', declared ' + p.declared_sex) : '')));
    }
    var r = role();
    if (r) { bits.push(el('span', {}, r)); }
    if (run.state === 'finished') {
      var hh = U.hours(run.started_utc, run.finished_utc);
      bits.push(el('span', {}, 'finished' + (run.finished_utc ? ' ' + U.dm(run.finished_utc) : '') + (hh ? ' after ' + hh.toFixed(1) + ' h' : '')));
    } else if (running()) {
      bits.push(el('span', {}, done + ' of ' + planned.length + ' steps done'));
    } else if (run.state === 'stopped') {
      bits.push(el('span', {}, 'stopped' + (run.last_heard_utc ? ', no news since ' + U.dm(run.last_heard_utc) : '') + ', ' + done + ' of ' + planned.length + ' steps done'));
    }
    bits.forEach(function (b) { add(line, b); });
    // Breadcrumb: people / family / person.
    var crumbs = document.getElementById('crumbs');
    var name = document.getElementById('run-name');
    var home = document.getElementById('home-crumb');
    crumbs.textContent = '';
    add(crumbs, home, el('span', { 'class': 'sep', 'aria-hidden': 'true' }, '/'));
    var fam = state && state.plugins && state.plugins.family;
    ((fam && fam.kind === 'member' && fam.families) || []).forEach(function (f) {
      add(crumbs, el('a', { href: f.url }, 'Family ' + f.name), el('span', { 'class': 'sep', 'aria-hidden': 'true' }, '/'));
    });
    add(crumbs, name);
    // The strip under the timeline: one cell per planned step.
    var strip = document.getElementById('strip');
    strip.textContent = '';
    planned.forEach(function (s) { add(strip, el('span', { 'class': 'cell st-' + s.state.replace(/\s+/g, '-'), title: s.step + ' ' + s.label + ': ' + s.state })); });
    var nrun = planned.filter(function (s) { return s.state === 'running'; }).length;
    var failed = planned.filter(function (s) { return s.state === 'failed'; }).length;
    var t = done + ' of ' + planned.length + ' steps done';
    if (run.state === 'finished') { t = 'Finished' + (run.finished_utc ? ' ' + U.day(run.finished_utc) : '') + ': ' + t; }
    else if (run.state === 'stopped') {
      t = 'Stopped' + (run.last_heard_utc ? ': no news since ' + U.day(run.last_heard_utc) : '') + ', ' + t +
        ', with no end record, so this report may be incomplete';
    } else if (nrun) { t = 'Running: ' + t + ', ' + nrun + ' running now'; }
    if (failed) { t += ', ' + failed + ' failed'; }
    var rc = state && state.plugins && state.plugins.recheck;
    if (rc && rc.checked_utc) { t += '. Databases re-checked ' + U.day(rc.checked_utc); }
    document.getElementById('strip-text').textContent = t + '.';
    runHero(planned, done, nrun, failed);
    root.setAttribute('data-pgp-run', run.state || '');
  }
  // While the run goes: the hero (31 of 54), the timing and what runs now.
  function runHero(planned, done, nrun, failed) {
    var hero = document.getElementById('runhero');
    var badge = document.getElementById('live-badge');
    hero.hidden = !running();
    if (badge) { badge.hidden = !running(); }
    var list = document.getElementById('running-list');
    var box = document.getElementById('running-now');
    list.textContent = '';
    var now = Date.parse((state && state.generated_utc) || '');
    // A stopped run's steps still read running, but nothing runs now.
    planned.filter(function (s) { return running() && s.state === 'running'; }).forEach(function (s) {
      var li = el('li', {});
      var el2 = isFinite(now) && s.started_utc ? (now - Date.parse(s.started_utc)) / 60e3 : null;
      var usual = s.hours ? s.hours * 60 : null;
      add(li, add(el('div', { 'class': 'rn' }), el('span', {}, s.label),
        el('span', {}, (el2 !== null ? Math.round(el2) + ' min' : '') + (usual ? ' of ~' + Math.round(usual) + ' min' : ''))));
      if (el2 !== null && usual) {
        var pb = el('div', { 'class': 'pb' });
        var i = el('i', {});
        i.style.setProperty('--p', Math.min(100, 100 * el2 / usual).toFixed(1) + '%');
        add(li, add(pb, i));
      }
      add(list, li);
    });
    box.hidden = !list.childNodes.length;
    if (!running()) { return; }
    var ops = state.plugins && state.plugins.ops;
    // hours_left is the ops plugin's {low_s, high_s} range from the usual step times, or null.
    var hl = ops && ops.hours_left;
    var left = hl && typeof hl.high_s === 'number' ? hl.high_s / 3600 : null;
    var low = hl && typeof hl.low_s === 'number' ? hl.low_s / 3600 : null;
    var span = left === null ? null : (low !== null && Math.round(low) !== Math.round(left) ? Math.max(0, Math.round(low)) + '–' + Math.max(1, Math.round(left)) + ' h' : '~' + Math.max(1, Math.round(left)) + ' h');
    var retries = ops && ops.memory ? ops.memory.retries || 0 : 0;
    var alerts = (ops && ops.alerts) || [];
    var lede = document.getElementById('run-lede');
    lede.textContent = '';
    add(lede, 'This report is being written as the steps finish, each section in its place. ',
      el('b', {}, done + ' of ' + planned.length + ' are done'),
      span !== null ? ', about ' + span.replace(' h', ' hours') + ' remain. ' : '. ',
      failed ? el('b', {}, failed + ' failed: the timeline says which. ') : '',
      alerts.length ? alerts.map(function (a) { return a.text; }).join('; ') + '. ' : (failed ? '' : 'Nothing needs you. '),
      findingsNow().length ? 'The first findings are already in.' : 'No section is written yet.');
    var big = document.getElementById('bignum');
    big.textContent = String(done);
    add(big, el('small', {}, 'of ' + planned.length + ' steps'));
    var bar = document.getElementById('stepbar');
    bar.textContent = '';
    [['d', done], ['r', nrun], ['f', failed]].forEach(function (x) {
      if (!x[1]) { return; }
      var i = el('i', { 'class': x[0] });
      i.style.setProperty('width', (100 * x[1] / Math.max(1, planned.length)).toFixed(2) + '%');
      add(bar, i);
    });
    var key = document.getElementById('stepkey');
    key.textContent = '';
    add(key, el('span', {}, done + ' done'), el('span', {}, nrun + ' running'),
      el('span', {}, (planned.length - done - nrun - failed) + ' waiting'), el('span', {}, failed + ' failed'));
    var timing = document.getElementById('timing');
    timing.textContent = '';
    var started = state.run.started_utc || (planned.map(function (s) { return s.started_utc; }).filter(Boolean).sort()[0]);
    var el3 = U.hours(started, state.generated_utc);
    [['Elapsed', el3 !== null ? el3.toFixed(1) + ' h' : 'not known', started ? 'since ' + U.hm(started) : ''],
      ['To go', span !== null ? span : 'not known', span !== null ? 'from usual step times' : 'no usual times yet'],
      ['Retries', String(retries), retries ? 'recovered unless a step failed' : 'none']].forEach(function (x) {
      add(timing, add(el('div', {}), el('span', { 'class': 'k' }, x[0]), el('span', { 'class': 'v' }, x[1]), el('span', { 'class': 'r' }, x[2])));
    });
  }

  // --- the route: #s/<card>, #s<n>, #f/<finding>, #q/<chapter>, #progress -----------------------------
  var CHAPTER_SECTION = { quality: 5, clinical: 8, carrier: 3, medicines: 2, ancestry: 6, traits: 7, other: 10, look: 1, leaves: 9 };
  // The element the address names, and whether it is a card; null when it names nothing on the page.
  function routeTarget() {
    var h;
    try { h = decodeURIComponent((window.location.hash || '').slice(1)); } catch (e) { return null; }
    if (!h) { return null; }
    var target = null, isCard = false;
    if (h.indexOf('s/') === 0 || (h && !/^(s\d+|f\/.*|q\/.*|progress)$/.test(h) && cardById(h))) {
      target = document.getElementById('card-' + (h.indexOf('s/') === 0 ? h.slice(2) : h));
      isCard = true;
    } else if (h.indexOf('f/') === 0) { var hf = hashFinding(); target = document.getElementById('f-' + (hf ? handle(hf) : h.slice(2))); }
    else if (h.indexOf('q/') === 0) { target = document.getElementById('s' + (CHAPTER_SECTION[h.slice(2)] || 1)); }
    else { target = document.getElementById(h); }
    return target ? { target: target, isCard: isCard } : null;
  }
  // Open the folds around the place the address names, and mark it when it is a card. No scroll.
  function openAround(t) {
    var target = t.target;
    if (target.id === 'progress') { var fold = document.getElementById('progress-fold'); if (fold) { fold.open = true; } }
    // Open every fold around the target (a card inside a folded section on a phone sits in two), and a
    // folded section's own body when the link names the section.
    for (var d = target.closest ? target.closest('details') : null; d; d = d.parentElement && d.parentElement.closest('details')) { d.open = true; }
    var own = target.matches && target.matches('section.sec') && target.querySelector(':scope > details.secfold');
    if (own) { own.open = true; }
    if (t.isCard) {
      Array.prototype.forEach.call(panel.querySelectorAll('.result.target'), function (x) { x.classList.remove('target'); });
      target.classList.add('target');
    }
  }
  function route() {
    syncHash();
    var t = routeTarget();
    if (!t) { return; }
    var target = t.target;
    openAround(t);
    target.scrollIntoView({ block: 'start' });
    if (target.focus) { target.focus({ preventScroll: true }); }
    markToc(target);
  }
  // The contents mark the section the address names; a redraw marks it again, so the mark never depends on
  // whether the page drew again after the address was read.
  function markToc(target) {
    Array.prototype.forEach.call(rail.querySelectorAll('[data-toc]'), function (a) {
      var sec = target.closest ? target.closest('section.sec') : null;
      a.setAttribute('aria-current', sec && sec.id === 's' + a.getAttribute('data-toc') ? 'true' : 'false');
    });
  }
  window.addEventListener('hashchange', route);
  // A link to the place the address already names fires no hashchange: open its folds all the same. A
  // click that opens a new tab or window (Ctrl, Cmd, Shift, Alt, a button other than the main one) leaves
  // this page as it is.
  document.addEventListener('click', function (e) {
    if (e.button !== 0 || e.ctrlKey || e.metaKey || e.shiftKey || e.altKey) { return; }
    var a = e.target && e.target.closest ? e.target.closest('a[href^="#"]') : null;
    if (a && a.getAttribute('href') === window.location.hash) { route(); }
  });
  // Turning a phone or resizing a window across 700 px draws the page again, folded or open.
  if (phoneQuery) {
    // The redraw keeps every fold as the reader left it, the phone's folds too when the page went wide
    // and back (foldState); open the ones around the place the address names as well, without scrolling
    // back to it (the reader may be elsewhere now).
    var redrawWide = function () { if (results) { redraw(true); var t = routeTarget(); if (t) { openAround(t); } } };
    if (phoneQuery.addEventListener) { phoneQuery.addEventListener('change', redrawWide); } else if (phoneQuery.addListener) { phoneQuery.addListener(redrawWide); }
  }

  // --- data ---------------------------------------------------------------------------------------
  function load() {
    var key = replay ? 'replay' : doneSet() + '|' + ((state && state.results_stamp) || '');
    if (key === fetchedFor) { return; }
    fetchedFor = key;
    U.getJSON('results.json').then(function (r) {
      if (key !== fetchedFor) { return; }
      if (!(r && (r.chapters || r.note))) { fetchedFor = null; return; }
      results = r;
      resultsKey = key;
      redraw(true);
      route();
    }).catch(function () { if (key === fetchedFor) { fetchedFor = null; } });
  }
  // Every fold of the report as the reader left it, open or closed, by a key that names the same fold in
  // the next drawing: its data-card (a finding's caveat "More" carries its finding's id), else the element
  // with an id around it, its class and its place among the folds of that class there (a card's "other
  // lines", its APOE warning, its table of detection limits). A fold the
  // current drawing does not have (a phone's fold while the page is wide) keeps its last state.
  var foldState = {};
  function foldKeys() {
    var seen = {};
    return Array.prototype.map.call(panel.querySelectorAll('details'), function (d) {
      var key = d.getAttribute('data-card');
      if (!key) {
        key = ((d.parentElement && d.parentElement.closest('[id]')) || {}).id + ' ' + d.className;
        seen[key] = (seen[key] || 0) + 1;
        key += ' ' + seen[key];
      }
      return [key, d];
    });
  }
  function keepFolds() { foldKeys().forEach(function (x) { foldState[x[0]] = x[1].open; }); }
  function restoreFolds() { foldKeys().forEach(function (x) { if (foldState[x[0]]) { x[1].open = true; } }); }
  function redraw(force) {
    var k = JSON.stringify([doneSet(), state && state.run && state.run.state, !!results, !!family, familyStatus, Object.keys(shown),
      state && state.plugins && state.plugins.recheck && state.plugins.recheck.mtime,
      // Every step's state: one can move from waiting to not finished with the done set unchanged.
      ((state && state.steps) || []).map(function (s) { return s.step + ':' + s.state; })]);
    if (!force && k === drawnKey) { header(); return; }
    drawnKey = k;
    var y = window.scrollY;
    keepFolds();
    header();
    drawToc();
    drawTabbar();
    drawLead();
    drawSections();
    restoreFolds();
    // A gated result's fold follows its switch, whatever fold state the reader left (setShown).
    Array.prototype.forEach.call(panel.querySelectorAll('details[data-gate]'), function (d) { d.open = !!shown[d.getAttribute('data-gate')]; });
    window.scrollTo(0, y);
    syncHash();
    var rt = routeTarget();
    if (rt) { markToc(rt.target); }
    document.getElementById('progress-q').textContent = running()
      ? 'Is it moving, when will it end, is anything wrong? One line per step.' : 'One line per step: when it ran, and what it found.';
    // The step timeline closes the report, folded once the page has a report to show, while the run goes
    // and after it ends: the page answers "is there anything serious?" first. The band under the header
    // (runHero) and the contents link to it. A page with no report yet keeps it open: it is all there is.
    var fold = document.getElementById('progress-fold');
    if (fold && !folded && sectionsList().length) { folded = true; fold.open = window.location.hash === '#progress'; }
  }

  window.PGP.register({
    name: 'person',
    onState: function (s) {
      state = s;
      loadFamily();
      if (!replay || !results) { load(); }
      redraw(false);
    }
  });
}());
