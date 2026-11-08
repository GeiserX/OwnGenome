// home.js: the live site's home page (assets/live/home.html), the Journal design: "does anyone have
// something to act on, and what is running?" From /runs.json (bin/owngenome_live.py runs_index): a
// headline sentence, the findings that need attention across everyone, the run being written first,
// then each family as a block (a small pedigree and one tile per person, with the person's lead and the
// counts by level), the people in no family, and the log of every run. Colour follows the person: each
// member of a family keeps one colour on this page. The plain list of runs is what shows without scripts.
(function () {
  'use strict';
  var U = window.PGPUI;
  if (!U) { return; }
  var el = U.el, add = U.add;
  var root = document.documentElement;
  var PCOL = ['var(--p1)', 'var(--p2)', 'var(--p3)', 'var(--p4)', 'var(--p5)'];
  var WORDS = ['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine'];
  function word(n) { return WORDS[n] || String(n); }
  function cap(t) { return t.charAt(0).toUpperCase() + t.slice(1); }

  function role(m, fam) {
    var kids = fam.members.filter(function (x) { return x.father === m.sample || x.mother === m.sample; });
    if (kids.length) { return { male: 'father', female: 'mother' }[m.sex] || 'parent'; }
    if (fam.members.some(function (x) { return x.sample === m.father || x.sample === m.mother; })) {
      return { male: 'son', female: 'daughter' }[m.sex] || 'child';
    }
    return m.partner ? 'partner' : 'member';
  }
  function progressText(r) {
    if (!r) { return 'no genome on this server'; }
    if (r.state === 'finished') { return 'finished' + (r.finished_utc ? ' ' + U.dm(r.finished_utc) : ''); }
    if (r.state === 'running') { return 'being written: ' + r.done + ' of ' + r.planned + ' steps'; }
    if (!r.planned) { return 'not started'; }
    return r.done + ' of ' + r.planned + ' steps';
  }
  // The counts by level a tile shows: never a bare "nothing": no count means no finding of that level
  // in what was assessed, and the tile says the report is the place to read it.
  function countsRow(r) {
    var row = el('div', { 'class': 'counts' });
    var lv = (r && r.levels) || null;
    if (!lv) { add(row, el('span', { 'class': 'dim' }, r && r.done ? 'Results so far on the page' : 'No result yet')); return row; }
    [['act', 'doctor'], ['medicine', 'prescriber'], ['review', 'look'], ['carrier', 'family'], ['info', 'info']].forEach(function (x) {
      if (lv[x[0]]) { add(row, U.level(x[0], lv[x[0]] + ' ' + x[1])); }
    });
    if (!row.childNodes.length) { add(row, el('span', { 'class': 'dim' }, 'No finding in what was assessed; read the caveats')); }
    return row;
  }
  function tile(m, fam, r, colour) {
    var div = el('div', { 'class': 'ptile person-card' + (r ? ' has-run' : ' norun') + (r && r.state === 'running' ? ' running' : ''),
      'data-sample': m.sample });
    if (colour) { div.style.setProperty('--pcol', colour); }
    add(div, el('span', { 'class': 'bar-top' }));
    var h = el('h3', {});
    add(h, r ? el('a', { href: r.url, 'class': 'person-link' }, m.sample) : el('span', {}, m.sample));
    add(div, h);
    var meta = el('div', { 'class': 'meta' });
    add(meta, (fam ? role(m, fam) + ' · ' : '') + progressText(r));
    if (r && r.failed) { add(meta, ' · ', el('span', { 'class': 'warn' }, r.failed + ' step' + (r.failed > 1 ? 's' : '') + ' not assessed')); }
    add(div, meta);
    var child = !!fam && fam.members.some(function (x) { return x.sample === m.father || x.sample === m.mother; });
    var lead = r && (child ? r.lead_closed || '' : r.lead);
    if (lead) {
      // Two sentences: the first often says only that nothing needs a doctor, the second names what differs.
      var one = U.firstSentence(lead), rest = lead.slice(one.length).trim();
      add(div, el('p', { 'class': 'lead' }, rest ? one + ' ' + U.firstSentence(rest) : one));
    }
    else if (!r) { add(div, el('p', { 'class': 'lead' }, 'In the family tree, with no genome on this server.')); }
    add(div, countsRow(r));
    if (r) { add(div, el('a', { 'class': 'read', href: r.url }, 'Read the report →')); }
    return div;
  }
  // A small pedigree, one row per generation; filled when sequenced, dashed when not. Square male,
  // circle female, diamond when the sex is not known.
  function miniPedigree(fam, colourOf) {
    var ms = fam.members, W = 260, GAP = 66;
    var by = {};
    ms.forEach(function (m) { by[m.sample] = m; });
    var depth = {};
    function gen(m, seen) {
      if (depth[m.sample] !== undefined) { return depth[m.sample]; }
      if (seen[m.sample]) { return 0; }  // a loop in a malformed file: stop rather than recurse
      seen[m.sample] = true;
      var ps = [m.father, m.mother].filter(function (p) { return by[p]; });
      depth[m.sample] = ps.length ? 1 + Math.max.apply(null, ps.map(function (p) { return gen(by[p], seen); })) : 0;
      return depth[m.sample];
    }
    ms.forEach(function (m) { gen(m, {}); });
    var rows = [];
    ms.forEach(function (m) { (rows[depth[m.sample]] = rows[depth[m.sample]] || []).push(m); });
    var H = 34 + (rows.length - 1) * GAP + 32;
    var s = U.svg('svg', { viewBox: '0 0 ' + W + ' ' + H, width: W, height: H, role: 'img', 'aria-label': 'Pedigree of ' + fam.name });
    var pos = {};
    rows.forEach(function (list, g) {
      (list || []).forEach(function (m, i) { pos[m.sample] = { x: (i + 1) * W / (list.length + 1), y: 34 + g * GAP }; });
    });
    ms.forEach(function (k) {
      var ps = [k.father, k.mother].filter(function (p) { return pos[p]; });
      if (!ps.length) { return; }
      var py = Math.max.apply(null, ps.map(function (p) { return pos[p].y; }));
      var mx = ps.reduce(function (a, p) { return a + pos[p].x; }, 0) / ps.length;
      var mid = (py + pos[k.sample].y) / 2;
      if (ps.length === 2) { s.appendChild(U.svg('line', { x1: pos[ps[0]].x, x2: pos[ps[1]].x, y1: py, y2: py, stroke: '#8f8d86' })); }
      s.appendChild(U.svg('line', { x1: mx, x2: mx, y1: py, y2: mid, stroke: '#8f8d86' }));
      s.appendChild(U.svg('line', { x1: mx, x2: pos[k.sample].x, y1: mid, y2: mid, stroke: '#8f8d86' }));
      s.appendChild(U.svg('line', { x1: pos[k.sample].x, x2: pos[k.sample].x, y1: mid, y2: pos[k.sample].y - 11, stroke: '#8f8d86' }));
    });
    ms.forEach(function (m) {
      var p = pos[m.sample], c = colourOf[m.sample] || '#8f8d86';
      var attrs = { fill: '#15161a', stroke: m.sequenced ? c : '#8f8d86', 'stroke-width': 2 };
      if (!m.sequenced) { attrs['stroke-dasharray'] = '3 3'; }
      var shape = m.sex === 'male' ? U.svg('rect', Object.assign({ x: p.x - 11, y: p.y - 11, width: 22, height: 22 }, attrs))
        : m.sex === 'female' ? U.svg('circle', Object.assign({ cx: p.x, cy: p.y, r: 11 }, attrs))
          : U.svg('polygon', Object.assign({ points: [p.x, p.y - 13, p.x + 13, p.y, p.x, p.y + 13, p.x - 13, p.y].join(' ') }, attrs));
      s.appendChild(shape);
      var t = U.svg('text', { x: p.x, y: p.y + 26, 'text-anchor': 'middle', fill: '#c3c2b7', 'font-size': 11 });
      t.textContent = m.sample;
      s.appendChild(t);
    });
    return s;
  }
  function rank(r) { return !r ? 3 : r.state === 'running' ? 0 : r.state === 'finished' ? 1 : 2; }
  function hoursOf(r) { var h = U.hours(r.started_utc, r.finished_utc); return h === null ? '—' : h.toFixed(1); }

  function draw(data) {
    var main = document.getElementById('home');
    main.textContent = '';
    var runs = data.runs || [], fams = data.families || [];
    var byName = {};
    runs.forEach(function (r) { byName[r.name] = r; });
    var colourOf = {};
    fams.forEach(function (f) { f.members.forEach(function (m, i) { if (!colourOf[m.sample]) { colourOf[m.sample] = PCOL[i % PCOL.length]; } }); });
    var going = runs.filter(function (r) { return r.state === 'running'; });
    var done = runs.filter(function (r) { return r.state === 'finished'; });
    var acts = runs.reduce(function (a, r) { return a + ((r.levels && r.levels.act) || 0); }, 0);
    var assessed = runs.filter(function (r) { return r.levels; }).length;
    // A child of a family: results to act on stay closed on this page too (the report opens them).
    var children = {};
    fams.forEach(function (f) {
      f.members.forEach(function (m) { if (f.members.some(function (x) { return x.sample === m.father || x.sample === m.mother; })) { children[m.sample] = true; } });
    });
    function closedTitle(r, f) { return children[r.name] && f.level === 'act' ? 'A result to act on, closed until opened on the report' : f.title; }
    // The dateline and the headline sentence.
    var dl = el('div', { 'class': 'dateline' });
    var today = new Date();
    add(dl, el('span', { 'class': 'd' }, today.toUTCString().slice(0, 16)),
      el('span', {}, runs.length + ' genome' + (runs.length === 1 ? '' : 's') + ' · ' + fams.length + ' famil' + (fams.length === 1 ? 'y' : 'ies') +
        ' · ' + going.length + ' report' + (going.length === 1 ? '' : 's') + ' being written'));
    add(main, dl);
    var hero = el('section', { 'class': 'homehero' });
    var left = el('div', {});
    add(left, el('div', { 'class': 'kicker' }, 'On this server'));
    var h1 = el('h1', { 'class': 'h-display' });
    if (!runs.length) { add(h1, 'No genome here yet.'); }
    else {
      if (going.length) { add(h1, cap(word(going.length)) + ' report' + (going.length > 1 ? 's are ' : ' is '), el('em', {}, 'being written'), '. '); }
      add(h1, (going.length ? cap(word(done.length)) + (done.length === 1 ? ' is' : ' are') : cap(word(done.length)) + ' report' + (done.length === 1 ? ' is' : 's are')) + ' finished, and ' +
        (acts ? word(acts) + ' finding' + (acts === 1 ? ' should' : 's should') + ' go to a doctor.'
          : !assessed ? 'none has results to read yet.'
            : assessed < runs.length ? 'no finding in the ' + word(assessed) + ' with results so far needs a doctor.'
              : 'no finding in what was assessed needs a doctor.'));
    }
    add(left, h1);
    var sub = fams.map(function (f) { return 'Family ' + f.name + ': ' + f.members.filter(function (m) { return m.sequenced; }).length + ' genomes, read alone and together.'; });
    if (going.length) { sub.push(going.map(function (r) { return r.name; }).join(', ') + ' ' + (going.length > 1 ? 'are' : 'is') + ' still running.'); }
    if (sub.length) { add(left, el('p', { 'class': 'sub' }, sub.join(' '))); }
    if (!runs.length) { add(left, el('p', { 'class': 'sub' }, 'Start one with scripts/run-all.sh and it shows here.')); }
    var right = el('div', {});
    add(right, el('div', { 'class': 'fig-label' }, 'Needs attention, across everyone'));
    var ul = el('ul', { 'class': 'attnlist', id: 'attention' });
    var order = { act: 0, medicine: 1, carrier: 2, review: 3, trust: -1 };
    var items = [];
    // One line per person and level: the first finding names the line, the rest are counted on it.
    runs.forEach(function (r) {
      var seen = {};
      (r.top || []).forEach(function (f) {
        if (order[f.level] === undefined) { return; }
        if (seen[f.level]) { return; }
        // The rest of the level from its count, not from the first four findings the index names.
        seen[f.level] = { r: r, f: f, more: Math.max(0, ((r.levels && r.levels[f.level]) || 1) - 1) };
        items.push(seen[f.level]);
      });
    });
    items.sort(function (a, b) { return order[a.f.level] - order[b.f.level]; });
    items.slice(0, 6).forEach(function (x) {
      var a = el('a', { href: x.r.url + 'finding/' + encodeURIComponent(x.f.id) + '/' });
      var who = el('span', { 'class': 'who' });
      var pc = el('span', { 'class': 'pc' }, x.r.name);
      if (colourOf[x.r.name]) { pc.style.setProperty('--pcol', colourOf[x.r.name]); }
      add(who, pc, ' · ' + U.LEVEL_WORDS[x.f.level].toLowerCase());
      add(a, U.findingLevel(x.f.level, ''), who, el('span', { 'class': 'what' }, closedTitle(x.r, x.f) + (x.more ? ' and ' + x.more + ' more' : '')), el('span', { 'class': 'arr', 'aria-hidden': 'true' }, '→'));
      add(ul, add(el('li', {}), a));
    });
    if (!items.length) { add(ul, add(el('li', {}), el('p', { 'class': 'soft' }, runs.length ? 'No finding to act on, for a prescriber or for family planning in what the finished steps assessed. Each report says what its steps do not see.' : 'Nothing yet.'))); }
    add(right, ul);
    add(hero, left, right);
    add(main, hero);
    // Being written: a run in progress comes first.
    going.forEach(function (r) {
      var head = el('div', { 'class': 'blockhead' });
      add(head, el('h2', {}, 'Being written'), el('span', { 'class': 'dim' }, 'a run in progress comes first'), el('a', { 'class': 'go', href: r.url }, 'Watch it →'));
      add(main, head);
      var card = el('div', { 'class': 'writing-card', 'data-run': r.name });
      var c1 = el('div', {});
      add(c1, el('div', { 'class': 'kicker' }, 'Being written'));
      add(c1, add(el('h3', {}), el('a', { href: r.url }, r.name)));
      add(c1, el('div', { 'class': 'meta dim' }, r.started_utc ? 'started ' + U.dm(r.started_utc) + ' ' + U.hm(r.started_utc) : ''));
      var big = el('div', { 'class': 'bignum' }, String(r.done));
      add(big, el('small', {}, 'of ' + r.planned + ' steps'));
      add(c1, big);
      var bar = el('div', { 'class': 'stepbar' });
      [['d', r.done], ['r', r.running], ['f', r.failed]].forEach(function (x) {
        if (!x[1]) { return; }
        var i = el('i', { 'class': x[0] });
        i.style.setProperty('width', (100 * x[1] / Math.max(1, r.planned)).toFixed(2) + '%');
        add(bar, i);
      });
      add(c1, bar, add(el('div', { 'class': 'stepkey' }), el('span', {}, r.done + ' done'), el('span', {}, r.running + ' running'),
        el('span', {}, (r.planned - r.done - r.running - r.failed) + ' waiting')));
      var c2 = el('div', {});
      add(c2, el('div', { 'class': 'fig-label' }, 'Found so far'));
      var fs = el('ul', { 'class': 'fs' });
      (r.top || []).forEach(function (f) { add(fs, add(el('li', {}), U.findingLevel(f.level, ''), ' ', closedTitle(r, f))); });
      if (!fs.childNodes.length) { add(fs, el('li', {}, r.lead ? U.firstSentence(r.lead) : 'No section is written yet: the first ones appear as their steps end.')); }
      add(c2, fs);
      var c3 = el('div', {});
      add(c3, el('div', { 'class': 'fig-label' }, 'Timing'));
      var el2 = U.hours(r.started_utc, new Date().toISOString());
      add(c3, add(el('div', { 'class': 'timing' }), add(el('div', {}), el('span', { 'class': 'k' }, 'Elapsed'),
        el('span', { 'class': 'v' }, el2 !== null ? el2.toFixed(1) + ' h' : '—'))));
      add(card, c1, c2, c3);
      add(main, card);
    });
    // The families.
    var inFamily = {};
    fams.forEach(function (f) {
      var head = el('div', { 'class': 'blockhead family-block family-head', 'data-family': f.name });
      var seq = f.members.filter(function (m) { return m.sequenced; }).length;
      add(head, add(el('h2', {}), el('a', { href: f.url }, 'Family ' + f.name)),
        el('span', { 'class': 'dim' }, seq + ' genome' + (seq === 1 ? '' : 's') + (f.members.length > seq ? ', ' + (f.members.length - seq) + ' without one' : '')),
        el('a', { 'class': 'go', href: f.url }, 'Family report →'));
      add(main, head);
      var rowx = el('div', { 'class': 'family-row' });
      var side = el('div', { 'class': 'fam-side' });
      var cols = {};
      f.members.forEach(function (m) { cols[m.sample] = { 'var(--p1)': '#d95926', 'var(--p2)': '#3987e5', 'var(--p3)': '#199e70', 'var(--p4)': '#c98500', 'var(--p5)': '#9085e9' }[colourOf[m.sample]]; });
      add(side, miniPedigree(f, cols));
      add(side, el('p', { 'class': 'soft' }, 'The family report reads the genomes together: who is related how, what the children inherited, and each couple’s carrier screen.'));
      add(side, el('a', { href: f.url }, 'Read the family report →'));
      var ppl = el('div', { 'class': 'fam-people' });
      f.members.slice().sort(function (a, b) { return rank(byName[a.sample]) - rank(byName[b.sample]); }).forEach(function (m) {
        inFamily[m.sample] = true;
        add(ppl, tile(m, f, byName[m.sample], colourOf[m.sample]));
      });
      add(rowx, side, ppl);
      add(main, rowx);
    });
    var loose = runs.filter(function (r) { return !inFamily[r.name]; }).sort(function (a, b) { return rank(a) - rank(b) || (a.name < b.name ? -1 : 1); });
    if (loose.length) {
      var lh = el('div', { 'class': 'blockhead family-block' });
      add(lh, el('h2', {}, fams.length ? 'Not in a family' : 'People'));
      add(main, lh);
      var lp = el('div', { 'class': 'fam-people' });
      loose.forEach(function (r) { add(lp, tile({ sample: r.name, sex: r.declared_sex }, null, r, null)); });
      add(main, lp);
    }
    // Every run: the log, for whoever runs the server.
    if (runs.length) {
      var eh = el('div', { 'class': 'blockhead' });
      add(eh, el('h2', {}, 'Every run'), el('span', { 'class': 'dim' }, 'the log, for whoever runs the server'));
      add(main, eh);
      var t = el('table', { 'class': 'runlog', id: 'runlog' });
      add(t, add(el('thead', {}), add(el('tr', {}), el('th', {}, 'Run'), el('th', {}, 'State'), el('th', { 'class': 'r' }, 'Steps'),
        el('th', { 'class': 'r' }, 'Hours'), el('th', {}, 'Finished'), el('th', {}, 'Notes'))));
      var tb = el('tbody', {});
      runs.slice().sort(function (a, b) { return rank(a) - rank(b) || (a.name < b.name ? -1 : 1); }).forEach(function (r) {
        var pc = el('span', { 'class': 'pc' }, r.name);
        if (colourOf[r.name]) { pc.style.setProperty('--pcol', colourOf[r.name]); }
        var stc = el('span', { 'class': 'st-' + r.state }, r.state === 'finished' ? '✓ finished' : r.state);
        var notes = r.failed ? r.failed + ' step' + (r.failed > 1 ? 's' : '') + ' failed: not assessed' : r.running ? r.running + ' running now' : '';
        add(tb, add(el('tr', {}), add(el('td', {}), add(el('a', { href: r.url }), pc)), add(el('td', {}), stc),
          el('td', { 'class': 'r' }, r.done + ' / ' + r.planned), el('td', { 'class': 'r' }, hoursOf(r)),
          el('td', {}, r.finished_utc ? U.dm(r.finished_utc) + ', ' + U.hm(r.finished_utc) : r.state === 'running' ? 'running' : '—'),
          el('td', {}, notes)));
      });
      add(t, tb);
      add(main, t);
    }
    root.setAttribute('data-pgp-people', String(main.querySelectorAll('.person-card').length));
    root.setAttribute('data-pgp-status', 'ready');
  }

  // While a run goes, the page asks again every 30 seconds, as the family page does.
  function start() {
    U.getJSON(root.getAttribute('data-data') || '/runs.json').then(function (d) {
      if (!d) { return; }
      draw(d);
      if ((d.runs || []).some(function (r) { return r.state === 'running'; })) { setTimeout(start, 30000); }
    }).catch(function () { root.setAttribute('data-pgp-status', 'offline'); setTimeout(start, 30000); });
  }
  if (document.readyState === 'loading') { document.addEventListener('DOMContentLoaded', start); } else { start(); }
}());
