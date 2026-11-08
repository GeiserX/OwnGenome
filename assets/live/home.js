// home.js: the live site's home page (assets/live/home.html): every person, grouped by family, from
// /runs.json (bin/owngenome_live.py runs_index). A card shows the pedigree symbol, the person's place
// in the family, the run's progress and how many cards the report flags; runs still going come first
// in each family. The plain list of runs in the page is what shows without scripts.
(function () {
  'use strict';
  var U = window.PGPUI;
  if (!U) { return; }
  var el = U.el, add = U.add;
  var root = document.documentElement;

  function role(m, fam) {
    var kids = fam.members.filter(function (x) { return x.father === m.sample || x.mother === m.sample; });
    if (kids.length) { return { male: 'father', female: 'mother' }[m.sex] || 'parent'; }
    if (fam.members.some(function (x) { return x.sample === m.father || x.sample === m.mother; })) {
      return { male: 'son', female: 'daughter' }[m.sex] || 'child';
    }
    return m.partner ? 'partner' : 'member';
  }
  function progressText(r) {
    if (!r) { return 'No run on this server'; }
    if (r.state === 'finished') { return 'Finished' + (r.finished_utc ? ' ' + U.day(r.finished_utc) : ''); }
    if (r.state === 'running') { return 'Running: ' + r.done + ' of ' + r.planned + ' steps'; }
    if (!r.planned) { return 'Not started'; }
    return r.done + ' of ' + r.planned + ' steps';
  }
  function lookText(r) {
    if (!r || !r.results) { return r && r.done ? 'Results so far on the page' : ''; }
    if (r.attention === null || r.attention === undefined) { return 'Results on the page'; }
    return r.attention ? r.attention + ' worth a look' : 'Nothing flagged; see the caveats';
  }
  function personCard(name, sex, sequenced, place, r) {
    var li = el('li', { 'class': 'person-card' + (r ? ' has-run' : '') + (r && r.state === 'running' ? ' running' : '') });
    var ratio = r && r.planned ? r.done / r.planned : null;
    var inner = r ? el('a', { href: r.url, 'class': 'person-link' }) : el('div', { 'class': 'person-link' });
    add(inner, U.avatar(sex, sequenced, ratio));
    var txt = el('div', { 'class': 'person-text' });
    add(txt, el('span', { 'class': 'person-name' }, name));
    if (place) { add(txt, el('span', { 'class': 'person-place' }, place)); }
    if (sequenced || r) {
      if (r && r.planned) {
        var bar = el('span', { 'class': 'bar', role: 'img', 'aria-label': r.done + ' of ' + r.planned + ' steps done' });
        var fill = el('span', { 'class': 'bar-fill' });
        fill.style.setProperty('--p', String(ratio));
        add(bar, fill);
        add(txt, bar);
      }
      add(txt, el('span', { 'class': 'person-state' }, progressText(r)));
      var lk = lookText(r);
      if (lk) { add(txt, el('span', { 'class': 'person-look' + (r && r.attention ? ' attn' : '') }, lk)); }
    } else {
      add(txt, el('span', { 'class': 'person-state' }, 'Not sequenced'));
    }
    add(inner, txt);
    add(li, inner);
    return li;
  }
  function rank(r) { return !r ? 3 : r.state === 'running' ? 0 : r.state === 'finished' ? 1 : 2; }

  function draw(data) {
    var main = document.getElementById('home');
    main.textContent = '';
    var byName = {};
    (data.runs || []).forEach(function (r) { byName[r.name] = r; });
    var inFamily = {};
    (data.families || []).forEach(function (f) {
      var block = el('section', { 'class': 'family-block', 'data-family': f.name });
      var h = el('header', { 'class': 'family-head' });
      var seq = f.members.filter(function (m) { return m.sequenced; }).length;
      var h2 = el('h2', {});
      add(h2, el('a', { href: f.url }, 'Family ' + f.name));
      add(h, h2, el('span', { 'class': 'dim' }, seq + ' of ' + f.members.length + ' sequenced'));
      add(h, el('a', { href: f.url, 'class': 'open' }, 'Open the family page'));
      add(block, h);
      var ul = el('ul', { 'class': 'people' });
      f.members.slice().sort(function (a, b) { return rank(byName[a.sample]) - rank(byName[b.sample]); }).forEach(function (m) {
        inFamily[m.sample] = true;
        add(ul, personCard(m.sample, m.sex, m.sequenced, role(m, f), byName[m.sample]));
      });
      add(block, ul);
      add(main, block);
    });
    var loose = (data.runs || []).filter(function (r) { return !inFamily[r.name]; })
      .sort(function (a, b) { return rank(a) - rank(b) || (a.name < b.name ? -1 : 1); });
    if (loose.length) {
      var block = el('section', { 'class': 'family-block' });
      add(block, add(el('header', { 'class': 'family-head' }), el('h2', {}, (data.families || []).length ? 'Not in a family' : 'People')));
      var ul = el('ul', { 'class': 'people' });
      loose.forEach(function (r) { add(ul, personCard(r.name, r.declared_sex, true, null, r)); });
      add(block, ul);
      add(main, block);
    }
    if (!loose.length && !(data.families || []).length) {
      add(main, el('p', { 'class': 'empty' }, 'No run under this folder yet. Start one with scripts/run-all.sh and it shows here.'));
    }
    root.setAttribute('data-pgp-people', String(main.querySelectorAll('.person-card').length));
    root.setAttribute('data-pgp-status', 'ready');
  }

  function start() {
    U.getJSON(root.getAttribute('data-data') || '/runs.json').then(function (d) { if (d) { draw(d); } })
      .catch(function () { root.setAttribute('data-pgp-status', 'offline'); });
  }
  if (document.readyState === 'loading') { document.addEventListener('DOMContentLoaded', start); } else { start(); }
}());
