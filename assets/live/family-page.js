// family-page.js: a family's page (assets/live/family.html), from family.json
// (bin/owngenome_live_family.py): the pedigree with each sequenced person a link to their page, the
// members with their run progress, what steps 38 and 39 planned, each couple's carrier screen in full
// (both carry, then X-linked, then one partner), who is related how, and each child's chromosome
// map. It draws what the steps wrote and computes nothing about genetics. Text through textContent,
// the drawing through createElementNS; images from this server only.
(function () {
  'use strict';
  var U = window.PGPUI;
  if (!U) { return; }
  var el = U.el, add = U.add, root = document.documentElement;
  var RESULT = { shared: 'Both carry', x_linked: 'X-linked', one_partner: 'One partner' };

  function section(title, cls) {
    var s = el('section', { 'class': 'card ' + (cls || '') });
    add(s, el('h2', {}, title));
    return s;
  }
  function memberUrl(f, s) { return f.runs && f.runs[s] ? f.runs[s].url : null; }

  // The pedigree drawn from the layout the server gives: a sequenced person with a run here is a link.
  function pedigree(f) {
    var L = f.layout, st = f.state;
    var box = el('figure', { 'class': 'pedigree-figure' });
    if (!L) { add(box, el('p', { 'class': 'dim' }, st.error || 'No pedigree.')); return box; }
    var s = U.svg('svg', { viewBox: '0 0 ' + L.width + ' ' + L.height, role: 'group', 'aria-label': 'Pedigree of ' + st.name,
      'class': 'pedigree-svg', preserveAspectRatio: 'xMidYMin meet' });
    var lines = U.svg('g', { 'class': 'ped-lines' });
    L.lines.forEach(function (l) { lines.appendChild(U.svg('line', { x1: l[0], y1: l[1], x2: l[2], y2: l[3] })); });
    s.appendChild(lines);
    s.style.maxWidth = Math.round(L.width * 1.1) + 'px';  // a two-person pedigree stays its own size
    var h = L.size / 2;
    L.people.forEach(function (p) {
      var url = memberUrl(f, p.sample);
      var g = U.svg(url ? 'a' : 'g', url ? { href: url, 'class': 'ped-person link', 'aria-label': p.sample + ': open their page' }
        : { 'class': 'ped-person', 'aria-label': p.sample + (p.sequenced ? ': no run on this server' : ': not sequenced') });
      var cls = 'ped-shape' + (p.sequenced ? ' filled' : '');
      var shape = p.sex === 'male' ? U.svg('rect', { x: p.x - h, y: p.y - h, width: L.size, height: L.size, 'class': cls })
        : p.sex === 'female' ? U.svg('circle', { cx: p.x, cy: p.y, r: h, 'class': cls })
          : U.svg('polygon', { points: [p.x, p.y - h, p.x + h, p.y, p.x, p.y + h, p.x - h, p.y].join(' '), 'class': cls });
      g.appendChild(shape);
      var t = U.svg('text', { x: p.x, y: p.y + h + 17, 'text-anchor': 'middle', 'class': 'ped-name' });
      t.textContent = p.sample;
      g.appendChild(t);
      var r = f.runs && f.runs[p.sample];
      if (r && r.attention) {
        var b = U.svg('text', { x: p.x + h + 4, y: p.y - h + 4, 'class': 'ped-flag' });
        b.textContent = '▲ ' + r.attention;
        g.appendChild(b);
      }
      var title = U.svg('title', {});
      title.textContent = p.sample + ', ' + (p.sex || 'sex not known') + (p.sequenced ? ', sequenced' : ', not sequenced');
      g.appendChild(title);
      s.appendChild(g);
    });
    add(box, s);
    add(box, el('figcaption', { 'class': 'dim' }, 'Square: male. Circle: female. Diamond: sex not known. Filled: sequenced. ' +
      'A person with a run on this server opens their page; ▲ counts the cards their report flags.'));
    var d = el('details', { 'class': 'pedigree-image' });
    add(d, el('summary', {}, 'The drawing as an image'));
    var img = el('img', { src: st.pedigree, alt: 'Pedigree of ' + st.name, 'class': 'pedigree' });
    add(d, img);
    add(box, d);
    return box;
  }

  function members(f) {
    var s = section('Members', 'members');
    var ul = el('ul', { 'class': 'people' });
    f.state.members.forEach(function (m) {
      var r = f.runs && f.runs[m.sample];
      var li = el('li', { 'class': 'person-card' + (r ? ' has-run' : '') });
      var inner = r ? el('a', { href: r.url, 'class': 'person-link' }) : el('div', { 'class': 'person-link' });
      add(inner, U.avatar(m.sex, m.sequenced, r && r.planned ? r.done / r.planned : null));
      var txt = el('div', { 'class': 'person-text' });
      add(txt, el('span', { 'class': 'person-name' }, m.sample));
      var rel = [];
      if (m.father || m.mother) { rel.push('child of ' + [m.father, m.mother].filter(Boolean).join(' and ')); }
      if (m.partner) { rel.push('partner of ' + m.partner); }
      if (rel.length) { add(txt, el('span', { 'class': 'person-place' }, rel.join('; '))); }
      var state = !m.sequenced ? 'Not sequenced' : !r ? 'No run on this server'
        : r.state === 'finished' ? 'Finished' : r.done + ' of ' + r.planned + ' steps';
      add(txt, el('span', { 'class': 'person-state' }, state));
      if (r && r.results) {
        add(txt, el('span', { 'class': 'person-look' + (r.attention ? ' attn' : '') },
          r.attention ? r.attention + ' worth a look' : 'Nothing flagged; see the caveats'));
      }
      add(inner, txt);
      add(li, inner);
      add(ul, li);
    });
    add(s, ul);
    return s;
  }

  function analyses(f) {
    var plan = f.state.plan || [];
    var p = el('p', { 'class': 'analyses' });
    if (!plan.length) { add(p, 'Analyses: none planned yet (step 38 writes family_plan.tsv).'); return p; }
    add(p, el('span', { 'class': 'label' }, 'Analyses: '));
    plan.forEach(function (r, i) {
      if (i) { add(p, ', '); }
      var word = r.decision === 'runs' ? r.analysis.replace(/_/g, ' ') : r.analysis.replace(/_/g, ' ') + ' (' + r.decision + ')';
      add(p, el('span', { 'class': 'plan-item d-' + r.decision, title: r.reason || '' }, word));
    });
    return p;
  }

  function couples(f) {
    var s = section('Couple carrier screen', 'couples');
    add(s, el('p', { 'class': 'caveat-lead' }, f.couple_note));
    if (!f.state.couples.length) { add(s, el('p', { 'class': 'dim' }, 'The pedigree names no partners.')); }
    f.state.couples.forEach(function (c) {
      var p = el('p', { 'data-couple': c.a + ' ' + c.b, 'class': 'couple-line' });
      if (!c.screened) {
        p.textContent = c.a + ' and ' + c.b + ': not screened yet (step 39).';
        add(s, p);
        return;
      }
      var rows = (f.couples || {})[c.a + ' ' + c.b] || [];
      var one = rows.filter(function (r) { return r.result === 'one_partner'; }).length;
      // A count of zero is not a clean result: the screen sees only what it lists as assessed.
      p.textContent = c.a + ' and ' + c.b + ': ' + c.shared + ' gene(s) where both carry a finding, ' + c.x_linked +
        ' X-linked, ' + one + ' where one partner carries a finding, unless the other partner carries something the ' +
        'screen does not see. ';
      if (c.report) { add(p, el('a', { href: c.report }, 'What the screen does not see')); }
      add(s, p);
      add(s, U.table({ columns: ['Result', 'Gene', 'Condition', 'Chance for a child', c.a, c.b, 'Why'],
        rows: rows.map(function (r) {
          var aFirst = r.partner_a === c.a;
          return [RESULT[r.result] || r.result, r.gene, r.condition, r.offspring_risk,
            aFirst ? r.findings_a : r.findings_b, aFirst ? r.findings_b : r.findings_a, r.explanation];
        }) }, 12));
    });
    return s;
  }

  function related(f) {
    var s = section('Who is related how', 'related');
    var rel = f.relationships || [], ibd = f.ibd_pairs || [];
    if (!rel.length && !ibd.length) { add(s, el('p', { 'class': 'dim' }, 'Not run yet (step 38).')); return s; }
    var byPair = {};
    ibd.forEach(function (r) { byPair[[r.sample_a, r.sample_b].sort().join(' ')] = r; });
    var rows = rel.map(function (r) {
      var i = byPair[[r.sample_a, r.sample_b].sort().join(' ')] || {};
      return [r.sample_a + ' and ' + r.sample_b, r.declared, r.relatedness, r.ibs0, r.ibs2, i.shared_cm || '.',
        i.degree || '.', r.verdict, r.reason];
    });
    ibd.forEach(function (i) {
      if (!rel.some(function (r) { return [r.sample_a, r.sample_b].sort().join(' ') === [i.sample_a, i.sample_b].sort().join(' '); })) {
        rows.push([i.sample_a + ' and ' + i.sample_b, i.declared, i.relatedness, '.', '.', i.shared_cm, i.degree, '.', '.']);
      }
    });
    add(s, U.table({ columns: ['Pair', 'Declared', 'Relatedness', 'IBS0', 'IBS2', 'Shared cM', 'Degree', 'Verdict', 'Reason'],
      rows: rows }, 20));
    if (!ibd.length) { add(s, el('p', { 'class': 'dim' }, 'Shared segments need step 26b’s panel: only the somalier half is here.')); }
    return s;
  }

  function children(f) {
    var s = section('Inherited stretches', 'children');
    if (!f.state.children.length) { add(s, el('p', { 'class': 'dim' }, 'The pedigree names no child.')); }
    f.state.children.forEach(function (c) {
      var d = el('details', { 'data-child': c.child });
      add(d, el('summary', {}, c.child + ': father ' + (c.father || 'not named') + ', mother ' + (c.mother || 'not named')));
      if (c.segments) {
        add(d, el('img', { src: c.segments, alt: 'Chromosome map of ' + c.child, 'class': 'segments' }));
      } else {
        add(d, el('p', { 'class': 'dim' }, 'No chromosome map yet: step 38 draws one for a child with a sequenced parent.'));
      }
      add(s, d);
    });
    return s;
  }

  function more(f) {
    var s = section('Plan and reports', 'plan');
    var plan = f.state.plan || [];
    if (plan.length) {
      add(s, U.table({ columns: ['Analysis', 'Members', 'Decision', 'Reason'],
        rows: plan.map(function (r) { return [r.analysis, r.members, r.decision, r.reason]; }) }, 20));
    }
    if (f.state.reports.length) {
      var ul = el('ul', { 'class': 'lines' });
      f.state.reports.forEach(function (r) { add(ul, add(el('li', {}), el('a', { href: r.url }, r.name))); });
      add(s, ul);
    }
    return s;
  }

  function draw(f) {
    var main = document.getElementById('family');
    main.textContent = '';
    var page = el('section', { 'class': 'family-page', 'data-plugin': 'family', 'aria-label': 'Family ' + f.state.name });
    if (f.state.error) {
      add(page, el('p', {}, f.state.error));
      add(main, page);
      root.setAttribute('data-pgp-status', 'ready');
      return;
    }
    var top = el('div', { 'class': 'family-top' });
    var ped = section('Pedigree', 'ped');
    add(ped, pedigree(f), analyses(f));
    add(top, ped, members(f));
    add(page, top);
    var grid = el('div', { 'class': 'family-grid' });
    add(grid, couples(f), related(f));
    add(page, grid, children(f));
    if ((f.state.plan || []).length || f.state.reports.length) { add(page, more(f)); }
    add(main, page);
    root.setAttribute('data-pgp-status', 'ready');
  }

  // The page asks again every 30 s and redraws only when family.json changed (a member's run moved on,
  // step 39 finished). A failed first request says so, with a retry; a later one keeps what is shown.
  var shown = null;
  function failed() {
    if (shown !== null) { return; }
    var main = document.getElementById('family');
    main.textContent = '';
    var b = el('button', { type: 'button', 'class': 'more' }, 'Try again');
    b.addEventListener('click', function () { main.textContent = ''; add(main, el('p', { 'class': 'dim' }, 'Loading the family.')); fetchFamily(); });
    add(main, el('p', {}, 'The family could not be loaded.'), b);
    root.setAttribute('data-pgp-status', 'offline');
  }
  function fetchFamily() {
    return U.getJSON('family.json').then(function (f) {
      if (!f) { failed(); return; }
      var k = JSON.stringify(f);
      if (k !== shown) { shown = k; draw(f); }
    }).catch(failed);
  }
  function start() {
    fetchFamily();
    window.setInterval(function () { if (!document.hidden) { fetchFamily(); } }, 30000);
  }
  if (document.readyState === 'loading') { document.addEventListener('DOMContentLoaded', start); } else { start(); }
}());
