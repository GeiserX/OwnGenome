// finding.js: one finding on a page of its own (assets/live/finding.html, the route
// /runs/<run>/finding/<id>/ of bin/owngenome_live_results.py). It answers "what do I do about it?":
// the explanation, the plan under it, then the evidence (the numbers, the frequency on a log strip, the
// reads with a link to the genome viewer), the caveat, the related findings and the family. Everything
// comes from the run's results.json (bin/live_results.py findings); nothing is computed here about
// genetics. A result to act on about a child of the family stays closed until it is opened on purpose.
(function () {
  'use strict';
  var U = window.PGPUI;
  if (!U) { return; }
  var el = U.el, add = U.add;
  var root = document.documentElement;
  var BASE = '../../';
  var fid = root.getAttribute('data-finding');
  var replay = root.getAttribute('data-mode') === 'replay';
  var main = document.getElementById('finding');
  var ROMAN = ['i', 'ii', 'iii', 'iv', 'v', 'vi', 'vii'];

  function abs(u) { return !u ? u : (u.charAt(0) === '/' ? u : BASE + u); }
  function note(head, body, cls) {
    var n = el('div', { 'class': 'note' + (cls ? ' ' + cls : '') });
    add(n, el('span', { 'class': 'nh' }, head), body);
    return n;
  }
  function dl(rows) {
    var d = el('dl', {});
    (rows || []).forEach(function (r) { add(d, el('dt', {}, r[0]), el('dd', {}, r[1])); });
    return d;
  }
  // The frequency on a log strip from 1 in a million to 1: where this change sits among people.
  function freqStrip(af) {
    var box = el('figure', {});
    add(box, el('div', { 'class': 'fig-label' }, 'How common it is (gnomAD), on a log scale'));
    var s = el('div', { 'class': 'fstrip', role: 'img', 'aria-label': af === null ? 'Frequency unknown' : 'Allele frequency ' + af });
    add(s, el('span', { 'class': 'axisline' }));
    [[1e-6, '1 in a million'], [1e-4, '1 in 10,000'], [1e-2, '1 in 100'], [1, 'everyone']].forEach(function (t) {
      var p = (Math.log10(t[0]) + 6) / 6 * 100;
      var k = el('span', { 'class': 'tick' }); k.style.setProperty('left', p + '%'); add(s, k);
      var l = el('span', { 'class': 'tl' }, t[1]); l.style.setProperty('left', p + '%'); add(s, l);
    });
    if (af !== null && af > 0) {
      var m = el('span', { 'class': 'me', title: 'This change: ' + af });
      m.style.setProperty('left', Math.max(0, Math.min(100, (Math.log10(af) + 6) / 6 * 100)) + '%');
      add(s, m);
    }
    add(box, s);
    add(box, el('figcaption', {}, af === null ? 'The frequency of this change is not known here, so it has not been checked: read the caveat.'
      : af === 0 ? 'Not seen in gnomAD.' : 'Allele frequency ' + af + '. A change above 1 in 100 is unlikely to be a high-penetrance finding.'));
    return box;
  }
  // A sketch of the reads at the position: what the genome viewer shows in full.
  function readsSketch(f) {
    var box = el('figure', { 'class': 'reads' });
    add(box, el('div', { 'class': 'fig-label' }, 'The reads at ' + f.viewer));
    var s = U.svg('svg', { viewBox: '0 0 680 90', role: 'img', 'aria-label': 'A sketch of reads over the position' });
    for (var i = 0; i < 9; i++) {
      var x = 30 + ((i * 97) % 260), w = 260 + ((i * 53) % 140);
      s.appendChild(U.svg('rect', { x: x, y: 6 + i * 9, width: w, height: 5, rx: 1, fill: 'rgba(255,255,255,.14)' }));
    }
    s.appendChild(U.svg('line', { x1: 340, x2: 340, y1: 0, y2: 90, stroke: '#c8ab78', 'stroke-dasharray': '3 3' }));
    box.appendChild(s);
    var cap = el('figcaption', {});
    add(cap, 'A sketch, not the data. ', el('a', { href: BASE + 'viewer/#' + encodeURIComponent(f.viewer) }, 'See the real reads in the genome viewer →'));
    add(box, cap);
    return box;
  }
  function draw(r, family) {
    var f = (r.findings || []).filter(function (x) { return x.id === fid; })[0];
    main.textContent = '';
    if (!f) {
      add(main, el('p', { 'class': 'empty' }, 'This finding is not in the results of this run now. The report lists every finding there is.'));
      root.setAttribute('data-pgp-status', 'ready');
      return;
    }
    // known, not family: the click that opens this finding draws again with family as it came, UNKNOWN kept.
    var unknown = family === UNKNOWN;
    var known = unknown ? null : family;
    var me = known && (known.state.members || []).filter(function (m) { return m.sample === r.sample; })[0];
    // A family that does not list this person cannot say they are not a child: closed as well.
    var child = unknown || (!!known && !me) || !!(me && (me.father || me.mother));
    // A finding that is not itself a result to act on but names one (a symptom ranking can rank the gene:
    // results.veil.findings, bin/live_results.py veil()) stays closed with it, as on the person page.
    // A tag is the list of every result it names: the finding stays closed while any of them is.
    var tags = (r.veil && r.veil.findings) || {};
    function tagOf(x) {
      if (x.level === 'act') { return null; }
      var t = tags[x.id];
      if (t === undefined) { t = ['*panel*']; }
      t = (Array.isArray(t) ? t : t ? [t] : []).filter(function (y) { return y !== x.id; });
      return t.length ? t : null;
    }
    // Opening this finding opens it and the results it names, nothing else: another result keeps its own veil.
    // A finding of a fixed-panel card (tag PANEL, bin/live_results.py FIXED_PANELS) names no result: it stays closed
    // while any result to act on is, whichever, as on the person page.
    var PANEL = '*panel*';
    var opened = {};
    if (draw.open) { opened[f.id] = true; (tagOf(f) || []).forEach(function (t) { if (t !== PANEL) { opened[t] = true; } }); }
    function actShut() { return (r.findings || []).some(function (x) { return x.level === 'act' && !opened[x.id]; }); }
    function shut(x) {
      if (!child) { return false; }
      if (x.level === 'act') { return !opened[x.id]; }
      return (tagOf(x) || []).some(function (t) { return t === PANEL ? !opened[x.id] && actShut() : !opened[t]; });
    }
    // Closed: no title, gene or variant anywhere on the page (its title, crumb, cover, rail) until opened.
    var closed = shut(f);
    var title = !closed ? f.title : f.level === 'act' ? 'A result to act on' : 'A finding about a result to act on';
    document.title = title + ': ' + r.sample + ': OwnGenome';
    document.getElementById('finding-crumb').textContent = title;
    var sec = (r.journal || []).filter(function (s) { return s.n === f.section; })[0];
    var same = (r.findings || []).filter(function (x) { return x.level === f.level; });
    var order = (r.findings || []).filter(function (x) { return x.level !== 'info'; });
    var next = order[order.indexOf(f) + 1];
    if (next && shut(next)) { next = null; }
    // The rail: back to the report, this page's parts, the next finding of the report.
    var rail = el('aside', { 'class': 'toc frail', 'aria-label': 'This finding' });
    var back = el('a', { 'class': 'back', href: BASE + '#f/' + f.id });
    add(back, el('span', { 'aria-hidden': 'true' }, '← '), r.sample + '’s report');
    add(rail, back);
    if (sec) { add(rail, el('div', { 'class': 'tsub' }, 'section ' + sec.n + ', ' + sec.title)); }
    add(rail, el('div', { 'class': 'group' }, 'This finding'));
    var parts = [['what', 'What it is'], ['plan-sec', 'What to do'], ['evidence', 'Evidence']];
    if (f.viewer) { parts.push(['reads', 'The reads']); }
    if (known) { parts.push(['family', 'In the family']); }
    parts.push(['related', 'Related results'], ['sources', 'Caveat and sources']);
    var tol = el('ol', {});
    parts.forEach(function (p, i) {
      var a = el('a', { href: '#' + p[0] });
      add(a, el('span', { 'class': 'tn' }, ROMAN[i]), el('span', {}, p[1]));
      add(tol, add(el('li', {}), a));
    });
    add(rail, tol);
    if (next) {
      add(rail, el('div', { 'class': 'group' }, 'Next in ' + r.sample + '’s report'));
      add(rail, el('a', { 'class': 'next', href: '../' + encodeURIComponent(next.id) + '/' }, U.LEVEL_WORDS[next.level] + ': ' + next.title + ' →'));
    }
    var cover = el('div', { 'class': 'fcover', id: 'what' });
    var lv = el('div', { 'class': 'fkick' });
    add(lv, U.findingLevel(f.level, f.level_words));
    add(lv, el('span', { 'class': 'dim' }, same.length === 1 ? 'The only finding at this level in ' + r.sample + '’s genome'
      : same.length + ' findings at this level in ' + r.sample + '’s genome'));
    add(cover, lv);
    add(cover, el('h1', { 'class': 'h-display' }, title));
    var chips = el('div', { 'class': 'chips' });
    var cs = [];
    if (f.gene) { cs.push(f.gene); }
    if (f.detail) { cs.push(f.detail); }
    if (f.viewer && cs.join(' ').indexOf(f.viewer) < 0) { cs.push(f.viewer); }
    cs.forEach(function (c) { add(chips, el('span', { 'class': 'chip mono' }, c)); });
    if (cs.length && !closed) { add(cover, chips); }
    var main0 = el('div', {});
    add(main0, cover);
    var acts = el('aside', { 'class': 'actions' });
    add(acts, el('a', { 'class': 'btn', href: BASE + '#f/' + f.id }, 'Back to the report'));
    add(acts, el('a', { 'class': 'btn ghost', href: BASE + 'report.html#' + encodeURIComponent(f.anchor || f.card) }, 'This section in the full report'));
    var body = el('div', { 'class': 'jlayout' });
    var col = el('div', {});
    var top = el('div', { 'class': 'row' });
    add(top, main0, acts);
    add(col, top);
    if (closed) {
      var v = el('div', { 'class': 'veil' });
      add(v, add(el('div', { 'class': 'vt' }), el('b', {}, 'A result about a child. '),
        (f.level === 'act' ? '' : 'This finding names a result to act on. ') +
        'A result to act on in adult life is usually not looked at before the child can choose. It stays closed until you open it.'));
      var b = el('button', { type: 'button', 'class': 'btn' }, 'Show the result');
      b.addEventListener('click', function () { draw.open = true; draw(r, family); });
      add(v, b);
      add(col, v);
      add(body, rail, col);
      add(main, body);
      root.setAttribute('data-pgp-status', 'ready');
      return;
    }
    var row1 = el('div', { 'class': 'row' });
    var m1 = el('div', {});
    add(m1, el('p', { 'class': 'lede' }, f.text));
    var box = el('div', { 'class': 'planbox lvl-' + f.level, id: 'plan-sec' });
    add(box, el('h2', { 'class': 'h3' }, 'What to do'));
    add(box, el('div', { 'class': 'secq' }, 'In this order. A research result: a certified laboratory confirms anything before it is acted on.'));
    var ol = el('ol', { 'class': 'todo', id: 'plan' });
    (f.plan || []).forEach(function (t) { add(ol, el('li', {}, t)); });
    add(box, ol);
    add(m1, box);
    var a1 = el('aside', {});
    if (f.level === 'carrier' && known) { add(a1, note('Family', add(el('span', {}), 'A carrier is not affected. ', el('a', { href: abs(familyUrl) }, 'The couple screen on the family page →')))); }
    add(row1, m1, a1);
    add(col, row1);
    // The evidence.
    var sec2 = el('section', { 'class': 'sec', id: 'evidence' });
    add(sec2, add(el('div', { 'class': 'sechead' }), el('span', { 'class': 'secnum' }, 'iii'), el('h2', { 'class': 'h2' }, 'The evidence')));
    add(sec2, el('div', { 'class': 'secq' }, 'The numbers behind it, as the run measured them.'));
    var row2 = el('div', { 'class': 'row' });
    var m2 = el('div', {});
    var facts = U.facts(f.evidence);
    if (facts) { add(m2, facts); } else { add(m2, el('p', { 'class': 'soft' }, 'The section in the full report holds the numbers.')); }
    if ('af' in f && (f.af !== null || f.card === 'clinvar')) { add(m2, freqStrip(typeof f.af === 'number' ? f.af : null)); }
    var a2 = el('aside', {});
    if (f.caveat && f.caveat.length) {
      var cv = el('div', {});
      f.caveat.forEach(function (t) { add(cv, el('p', {}, t)); });
      add(a2, note('Caveat', cv, 'caveat-note'));
    }
    add(row2, m2, a2);
    add(sec2, row2);
    add(col, sec2);
    if (f.viewer) {
      var sr = el('section', { 'class': 'sec', id: 'reads' });
      add(sr, add(el('div', { 'class': 'sechead' }), el('span', { 'class': 'secnum' }, 'iv'), el('h2', { 'class': 'h2' }, 'The reads')));
      add(sr, el('div', { 'class': 'secq' }, 'The sequenced fragments over the position, in the genome viewer.'));
      add(sr, add(el('div', { 'class': 'row' }), add(el('div', {}), readsSketch(f)), el('aside', {})));
      add(col, sr);
    }
    // Related findings and the family.
    var rel = (f.related || []).map(function (id) { return (r.findings || []).filter(function (x) { return x.id === id; })[0]; })
      .filter(function (x) { return x && !shut(x); });
    if (known) {
      var sf = el('section', { 'class': 'sec', id: 'family' });
      add(sf, add(el('div', { 'class': 'sechead' }), el('span', { 'class': 'secnum' }, ROMAN[parts.length - 3]), el('h2', { 'class': 'h2' }, 'In the family')));
      var fam = el('p', { 'class': 'soft' });
      add(fam, 'Who else in the family has a genome here: ');
      (known.state.members || []).filter(function (m) { return m.sample !== r.sample; }).forEach(function (m, i) {
        if (i) { add(fam, ', '); }
        var u = known.runs && known.runs[m.sample] && known.runs[m.sample].url;
        add(fam, u ? el('a', { href: abs(u) }, m.sample) : el('span', {}, m.sample + ' (no genome)'));
      });
      add(fam, '. ', el('a', { href: abs(familyUrl) }, 'The family report →'));
      add(sf, fam);
      add(col, sf);
    }
    var sec3 = el('section', { 'class': 'sec', id: 'related' });
    add(sec3, add(el('div', { 'class': 'sechead' }), el('span', { 'class': 'secnum' }, ROMAN[parts.length - 2]), el('h2', { 'class': 'h2' }, 'Related results')));
    add(sec3, el('div', { 'class': 'secq' }, 'Other findings of the same gene or the same step, and the family.'));
    var ul = el('ul', { 'class': 'list-findings' });
    rel.forEach(function (x) {
      var li = el('li', {});
      add(li, el('a', { 'class': 't', href: '../' + encodeURIComponent(x.id) + '/' }, x.title), U.findingLevel(x.level, x.level_words));
      add(li, el('span', { 'class': 'x' }, x.text));
      add(ul, li);
    });
    if (!rel.length) { add(ul, add(el('li', {}), el('span', { 'class': 'x' }, 'No other finding of this gene or this step.'))); }
    add(sec3, ul);
    add(col, sec3);
    var sec4 = el('section', { 'class': 'sec', id: 'sources' });
    add(sec4, add(el('div', { 'class': 'sechead' }), el('span', { 'class': 'secnum' }, ROMAN[parts.length - 1]), el('h2', { 'class': 'h2' }, 'Caveat and sources')));
    var src = el('p', { 'class': 'soft' });
    add(src, 'From the run’s ' + f.card + ' section (step ' + (f.steps || []).join(', ') + '). ',
      el('a', { href: BASE + 'report.html#' + encodeURIComponent(f.anchor || f.card) }, 'Open it in the full report'), '.');
    add(sec4, src);
    add(col, sec4);
    add(body, rail, col);
    add(main, body);
    root.setAttribute('data-pgp-status', 'ready');
  }
  var familyUrl = null;
  var UNKNOWN = { unknown: true };
  function familyOf(st) {
    // No state (state.json or frames.json refused, missing or not JSON), a state with no plugins, or a family
    // plugin that found a family report it could not read: the roles are unknown, never a standalone adult.
    if (!st || !st.plugins) { return Promise.resolve(UNKNOWN); }
    var fam = st.plugins.family;
    if (fam && fam.unknown) { return Promise.resolve(UNKNOWN); }
    var f = fam && fam.kind === 'member' && (fam.families || [])[0];
    if (!f) { return Promise.resolve(null); }
    familyUrl = f.url;
    // A family that does not answer: roles unknown, so a result to act on stays closed (draw()).
    return U.getJSON(abs(f.url) + 'family.json').then(function (d) { return d && d.state ? d : UNKNOWN; }).catch(function () { return UNKNOWN; });
  }
  function failed() {
    main.textContent = '';
    add(main, el('p', { 'class': 'empty' }, 'The results of this run could not be read just now.'));
    var b = el('button', { type: 'button', 'class': 'btn' }, 'Try again');
    b.addEventListener('click', function () { main.textContent = ''; add(main, el('p', { 'class': 'dim' }, 'Loading the finding.')); start(); });
    add(main, b);
    root.setAttribute('data-pgp-status', 'offline');
  }
  function start() {
    var st = replay ? U.getJSON(BASE + 'frames.json').then(function (d) { var fr = (d && d.frames) || []; return fr.length ? fr[fr.length - 1].state : null; })
      : U.getJSON(BASE + 'state.json');
    Promise.all([U.getJSON(BASE + 'results.json'), st.then(familyOf).catch(function () { return UNKNOWN; })]).then(function (x) {
      if (x[0]) { draw(x[0], x[1]); } else { failed(); }
    }).catch(failed);
  }
  if (document.readyState === 'loading') { document.addEventListener('DOMContentLoaded', start); } else { start(); }
}());
