// family-page.js: a family's page (assets/live/family.html), from family.json
// (bin/owngenome_live_family.py), the Journal design: a report in chapters, with a lead that answers
// "what does it mean for our children?". The pedigree with each sequenced person a link to their page, the
// members with their run progress, what steps 38 and 39 planned, each couple's carrier screen in full
// (both carry, then X-linked, then one partner), who is related how as a grid (a list of pairs on a
// phone), each child's stretches painted by parent, the de novo and compound-het candidates, and the
// members side by side. It draws what the steps wrote and computes nothing about genetics. Text through textContent,
// the drawing through createElementNS; images from this server only.
// A child's results to act on stay closed here as on the child's own page (person.js): each child's results.json
// is read by the shared gate, which passes no report data of a child with an unopened act to this view.
(function () {
  'use strict';
  window.PGPUI.ready.then(function () {
  if (window.PGPUI.gateOwnsPage()) { return; }
  var U = window.PGPUI;
  if (!U) { return; }
  var el = U.el, add = U.add, root = document.documentElement;
  var RESULT = { shared: 'Both carry', x_linked: 'X-linked', one_partner: 'One partner' };
  var STALE = ' (older than its input)';

  // A chapter: its roman number in brass, its title in the serif.
  var CHAPTER = { 'Pedigree': 'I', 'Who is related how': 'II', 'Couple carrier screen': 'III', 'Each child': 'IV',
    'Side by side': 'V', 'Plan and reports': 'App.' };
  function section(title, cls) {
    var s = el('section', { 'class': 'card ' + (cls || ''), id: 'ch-' + (cls || 'x') });
    var h = el('h2', {});
    if (CHAPTER[title]) { add(h, el('span', { 'class': 'secnum' }, CHAPTER[title])); }
    add(h, title);
    add(s, h);
    return s;
  }
  // Colour follows the person: the members in their order, one colour each, on every family figure.
  var PCOL = ['#d95926', '#3987e5', '#199e70', '#c98500', '#9085e9'];
  function colourOf(f, sample) {
    var i = (f.state.members || []).map(function (m) { return m.sample; }).indexOf(sample);
    return i >= 0 ? PCOL[i % PCOL.length] : null;
  }
  // The four equally likely children of two carriers of one recessive gene: 1 affected, 2 carriers, 1 neither.
  // Drawn only when the screen's own chance is 1 in 4: another chance (a partner with two copies) has
  // another split, and the row's words say it.
  function fourChildren(risk) {
    if (!/^\s*1 in 4\b/.test(risk || '')) { return null; }
    var box = el('div', { 'class': 'kids4' });
    var g = el('div', { 'class': 'grid4', role: 'img', 'aria-label': 'Each pregnancy: 1 in 4 affected, 2 in 4 carriers, 1 in 4 neither' });
    add(g, el('span', { 'class': 'aff', title: 'affected' }), el('span', { title: 'carrier' }), el('span', { title: 'carrier' }),
      el('span', { 'class': 'non', title: 'neither' }));
    add(box, g, add(el('div', { 'class': 'kt' }), el('b', {}, '1 in 4'), ' per pregnancy affected', el('br', {}),
      '2 in 4 carriers, not affected; 1 in 4 neither'));
    return box;
  }
  // The family's lead: what the screen, the relationship check and the children say, in plain words.
  function lead(f) {
    var st = f.state, out = [];
    var seq = (st.members || []).filter(function (m) { return m.sequenced; }).length;
    out.push([seq + ' genome' + (seq === 1 ? '' : 's') + ' of this family are read together here. ', false]);
    (st.couples || []).forEach(function (c) {
      if (!c.screened) { out.push([c.a + ' and ' + c.b + ' are not screened as a couple yet. ', false]); return; }
      if (c.shared) { out.push([c.a + ' and ' + c.b + ' both carry a change in ' + c.shared + ' gene' + (c.shared > 1 ? 's' : '') + ': chapter III gives the chance per pregnancy for each. ', true]); }
      else { out.push([c.a + ' and ' + c.b + ' share no carrier gene in what the screen sees; ', false], ['it does not see everything. ', false]); }
      if (c.x_linked) { out.push([c.x_linked + ' X-linked finding' + (c.x_linked > 1 ? 's' : '') + ' bear on sons. ', true]); }
    });
    var flags = (f.relationships || []).filter(function (r) { return r.verdict === 'mismatch' || r.verdict === 'warn'; }).length;
    var rels = f.relationships || [];
    var checked = rels.filter(function (r) { return r.verdict && r.verdict !== 'not_checked'; }).length;
    if (rels.length) {
      out.push([flags ? 'A measured relationship differs from the family tree: read chapter II before anything else. '
        : checked === rels.length ? 'The measured relationships match the family tree. '
          : 'Some relationships could not be checked: chapter II says which. ', !!flags]);
    }
    if ((st.children || []).length) {
      out.push([(st.children.length === 1 ? 'One child' : st.children.length + ' children') + ': chapter IV shows what each inherited from whom.', false]);
    }
    return out;
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
      var pc = p.sequenced ? colourOf(f, p.sample) : null;
      if (pc) { shape.style.setProperty('stroke', pc); }
      if (!p.sequenced) { shape.style.setProperty('stroke-dasharray', '4 3'); }
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
      var state = r && r.gated ? U.closedText : !m.sequenced ? 'Not sequenced' : !r ? 'No run on this server'
        : r.state === 'finished' ? 'Finished' : (r.state === 'stopped' ? 'Stopped, ' : '') + r.done + ' of ' + r.planned + ' steps';
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
      rows.filter(function (r) { return r.result === 'shared'; }).forEach(function (r) {
        var box = el('div', { 'class': 'mini', 'data-shared': r.gene });
        add(box, el('h3', { 'class': 'h3' }, (r.condition || r.gene) + ' (' + r.gene + ')'));
        var four = fourChildren(r.offspring_risk);
        // A chance other than 1 in 4 means a partner with two changed copies: not the carrier pattern.
        add(box, el('p', { 'class': 'soft' }, four ? 'Both carry a change in ' + r.gene + '. Carriers are not affected; a child who inherits both changes is.'
          : 'Both carry a change in ' + r.gene + ', and at least one of them has two changed copies, which is not the carrier pattern. A child who inherits a changed copy from each parent is affected.'));
        add(box, four || el('p', {}, 'Chance per pregnancy, as the screen gives it: ' + (r.offspring_risk || 'not given') + '.'));
        add(s, box);
      });
      add(s, U.table({ columns: ['Result', 'Gene', 'Condition', 'Chance for a child', c.a, c.b, 'Why'],
        rows: rows.map(function (r) {
          var aFirst = r.partner_a === c.a;
          return [RESULT[r.result] || r.result, r.gene, r.condition, r.offspring_risk,
            aFirst ? r.findings_a : r.findings_b, aFirst ? r.findings_b : r.findings_a, r.explanation];
        }) }, 12));
    });
    return s;
  }

  function pairKey(a, b) { return [a, b].sort().join(' '); }
  // Step 38's verdicts: ok (shown as matching), warn, mismatch, not_checked (too few sites to tell).
  function mismatch(v) { return v === 'mismatch' || v === 'warn'; }
  function related(f) {
    var s = section('Who is related how', 'related');
    var rel = f.relationships || [], ibd = f.ibd_pairs || [];
    if (!rel.length && !ibd.length) { add(s, el('p', { 'class': 'dim' }, 'Not run yet (step 38).')); return s; }
    var byRel = {}, byIbd = {}, names = [];
    function seen(n) { if (names.indexOf(n) < 0) { names.push(n); } }
    rel.forEach(function (r) { byRel[pairKey(r.sample_a, r.sample_b)] = r; seen(r.sample_a); seen(r.sample_b); });
    ibd.forEach(function (r) { byIbd[pairKey(r.sample_a, r.sample_b)] = r; seen(r.sample_a); seen(r.sample_b); });
    var order = (f.state.members || []).map(function (m) { return m.sample; });
    names.sort(function (a, b) { return order.indexOf(a) - order.indexOf(b); });
    function cellText(a, b) {
      var r = byRel[pairKey(a, b)] || {}, i = byIbd[pairKey(a, b)] || {};
      var lines = [];
      lines.push(i.degree || (r.relatedness ? 'relatedness ' + r.relatedness : 'not measured'));
      if (i.shared_cm) { lines.push(i.shared_cm + ' cM shared'); }
      var declared = r.declared || i.declared;
      if (declared) { lines.push('declared ' + declared); }
      if (r.verdict === 'not_checked') { lines.push('not checked: too few sites'); }
      else if (!r.verdict) { lines.push('relationship check not run'); }
      return { lines: lines, flag: mismatch(r.verdict), verdict: r.verdict, reason: r.reason };
    }
    // The grid: one cell per pair, the degree first, the declared relation and a flag when they differ.
    var grid = el('div', { 'class': 'rel-grid-wrap' });
    var t = el('table', { 'class': 'rel-grid' });
    var head = el('tr', {});
    add(head, el('th', { scope: 'col' }, ''));
    names.forEach(function (n) { add(head, el('th', { scope: 'col' }, n)); });
    add(t, add(el('thead', {}), head));
    var body = el('tbody', {});
    names.forEach(function (a) {
      var tr = el('tr', {});
      add(tr, el('th', { scope: 'row' }, a));
      names.forEach(function (b) {
        if (a === b) { add(tr, el('td', { 'class': 'self', 'aria-label': 'the same person' }, '·')); return; }
        if (!byRel[pairKey(a, b)] && !byIbd[pairKey(a, b)]) { add(tr, el('td', { 'class': 'none' }, 'not measured')); return; }
        var c = cellText(a, b);
        var td = el('td', { 'class': c.flag ? 'flag' : '', title: c.reason || '' });
        c.lines.forEach(function (l, i) { add(td, el('span', { 'class': i ? 'sub' : 'main' }, l)); });
        if (c.flag) { add(td, el('span', { 'class': 'chip chip-attn' }, '⚑ ' + c.verdict)); }
        add(tr, td);
      });
      add(body, tr);
    });
    add(t, body);
    add(grid, t);
    add(s, grid);
    // On a phone: a list of pairs.
    var ul = el('ul', { 'class': 'rel-pairs' });
    names.forEach(function (a, i) {
      names.slice(i + 1).forEach(function (b) {
        if (!byRel[pairKey(a, b)] && !byIbd[pairKey(a, b)]) { return; }
        var c = cellText(a, b);
        var li = el('li', { 'class': c.flag ? 'flag' : '' });
        add(li, el('span', { 'class': 'pair' }, a + ' and ' + b + ': '), c.lines.join(', '));
        if (c.flag) { add(li, ' ', el('span', { 'class': 'chip chip-attn' }, '⚑ ' + c.verdict + (c.reason ? ': ' + c.reason : ''))); }
        add(ul, li);
      });
    });
    add(s, ul);
    add(s, el('p', { 'class': 'dim legend-line' }, 'Each cell: the degree the shared stretches point to (or somalier’s ' +
      'relatedness), the centimorgans shared and the declared relation. ⚑ marks a pair the measurement does not match.'));
    var rows = rel.map(function (r) {
      var i = byIbd[pairKey(r.sample_a, r.sample_b)] || {};
      return [r.sample_a + ' and ' + r.sample_b, r.declared, r.relatedness, r.ibs0, r.ibs2, i.shared_cm || '.',
        i.degree || '.', r.verdict, r.reason];
    });
    ibd.forEach(function (i) {
      if (!byRel[pairKey(i.sample_a, i.sample_b)]) {
        rows.push([i.sample_a + ' and ' + i.sample_b, i.declared, i.relatedness, '.', '.', i.shared_cm, i.degree, '.', '.']);
      }
    });
    var d = el('details', { 'class': 'w-fallback' });
    add(d, el('summary', {}, 'Every measured number'));
    add(d, U.table({ columns: ['Pair', 'Declared', 'Relatedness', 'IBS0', 'IBS2', 'Shared cM', 'Degree', 'Verdict', 'Reason'],
      rows: rows }, 20));
    add(s, d);
    if (!ibd.length) { add(s, el('p', { 'class': 'dim' }, 'Shared segments need step 26b’s panel: only the somalier half is here.')); }
    return s;
  }

  // A child's stretches, painted: per chromosome the father's copy above the mother's, each stretch
  // coloured by the grandparent it came from (or by the parent when no grandparent was sequenced).
  var SIDE = { paternal: 'father', maternal: 'mother' };
  function painting(f, child, all) {
    var rows = all;
    var lengths = f.chrom_lengths || [];
    var longest = Math.max.apply(null, lengths.map(function (x) { return x[1]; }).concat([1]));
    var colours = Object.create(null), next = 1;
    function colourOf(r) {
      var who = r.ancestor && r.ancestor !== '.' ? r.ancestor : (r.parent && r.parent !== '.' ? r.parent : r.side + ' side');
      if (!colours[who]) { colours[who] = { n: next > 5 ? 5 : next++, side: r.side, basis: r.basis }; }
      return { who: who, n: colours[who].n };
    }
    var fig = el('figure', { 'class': 'painting', 'data-painting': child });
    var grid = el('div', { 'class': 'paint-grid', role: 'img', 'aria-label': 'Stretches ' + child + ' got from each parent, ' +
      'per chromosome; the table below lists them' });
    lengths.forEach(function (cl) {
      var chrom = cl[0];
      add(grid, el('span', { 'class': 'paint-chr' }, chrom.replace(/^chr/, '')));
      var both = el('span', { 'class': 'paint-tracks' });
      ['paternal', 'maternal'].forEach(function (side) {
        var track = el('span', { 'class': 'paint-track ' + side });
        track.style.setProperty('--len', (100 * cl[1] / longest).toFixed(2) + '%');
        rows.filter(function (r) { return r.chrom === chrom && r.side === side; }).forEach(function (r) {
          var a = parseInt(r.start, 10), b = parseInt(r.end, 10);
          if (!isFinite(a) || !isFinite(b) || b < a) { return; }
          var c = colourOf(r);
          var seg = el('span', { 'class': 'paint-seg anc' + c.n + (r.basis === 'other copy' ? ' other' : ''),
            title: SIDE[side] + '’s copy, ' + chrom + ':' + a + '-' + b + ', from ' + c.who + ' (' + r.sites + ' sites, ' +
              r.discordant + ' discordant)' });
          seg.style.setProperty('--a', (100 * a / cl[1]).toFixed(3) + '%');
          seg.style.setProperty('--w', Math.max(0.4, 100 * (b - a) / cl[1]).toFixed(3) + '%');
          add(track, seg);
        });
        add(both, track);
      });
      add(grid, both);
    });
    add(fig, grid);
    var leg = el('figcaption', { 'class': 'paint-legend' });
    add(leg, el('span', { 'class': 'dim' }, 'Upper line: the father’s copy; lower: the mother’s. Colour: who the stretch came from. '));
    Object.keys(colours).forEach(function (who) {
      var it = el('span', { 'class': 'paint-key' });
      add(it, el('span', { 'class': 'paint-swatch anc' + colours[who].n + (colours[who].basis === 'other copy' ? ' other' : '') }),
        who + (colours[who].basis === 'other copy' ? ' (the copy not from the sequenced parent)' : ''));
      add(leg, it);
    });
    add(leg, el('span', { 'class': 'dim' }, ' Grey: no stretch called (no informative site, a gap of the reference).'));
    add(fig, leg);
    var d = el('details', { 'class': 'w-fallback' });
    add(d, el('summary', {}, 'Every stretch (' + rows.length + ')'));
    add(d, U.table({ columns: ['Copy', 'Chromosome', 'Start', 'End', 'From', 'Sites', 'Discordant', 'Basis'],
      rows: rows.map(function (r) { return [SIDE[r.side] || r.side, r.chrom, r.start, r.end,
        r.ancestor !== '.' ? r.ancestor : r.parent, r.sites, r.discordant, r.basis]; }) }, 12));
    add(fig, d);
    return fig;
  }

  function children(f) {
    var s = section('Each child', 'children');
    if (!f.state.children.length) { add(s, el('p', { 'class': 'dim' }, 'The pedigree names no child.')); }
    // Whether step 38 assessed this child (family.json trio_status): "None listed" only for a child it assessed;
    // otherwise why not, and a file that links out of the family folder is named as one, never followed.
    function empty(child, key, file) {
      var st = ((f.trio_status || {})[child] || {})[key];
      if (st && st.assessed) { return 'Assessed (' + st.why + '): none listed in ' + file + ' for this child.'; }
      if (st && st.assessed === null) { return 'Not known whether assessed for this child: ' + st.why + '.'; }
      var why = st ? st.why : file + ' is not here';
      return 'Not assessed for this child: ' + why + '.';
    }
    f.state.children.forEach(function (c, i) {
      if (c.gated) {
        add(s, el('h3', {}, c.child), U.familyControls(c.child, fetchFamily));
        return;
      }
      var d = el('details', { 'data-child': c.child });
      if (i === 0) { d.open = true; }
      add(d, el('summary', {}, c.child + ': father ' + (c.father || 'not named') + ', mother ' + (c.mother || 'not named')));
      add(d, el('h3', {}, 'Inherited stretches'));
      var rows = (f.segments || {})[c.child];
      if (rows && rows.length) {
        add(d, painting(f, c.child, rows));
      } else if (c.segments) {
        // The drawing as an image is coloured by grandparent too, and cannot close one stretch: it stays closed
        // whole, as the painting does.
        add(d, el('img', { src: c.segments, alt: 'Chromosome map of ' + c.child, 'class': 'segments' }));
      } else {
        add(d, el('p', { 'class': 'dim' }, 'No chromosome map yet: step 38 draws one for a child with a sequenced parent.'));
      }
      var dn = (f.de_novo || []).filter(function (r) { return r.child === c.child; });
      add(d, el('h3', {}, 'De novo candidates'));
      if (dn.length) {
        add(d, el('p', { 'class': 'dim' }, 'In the child and in neither parent. Most candidates in a genome are genotype errors: ' +
          'check the reads before reading anything into one.'));
        // Every row closes while a result is: a row closed on its own would say where the result lies.
          add(d, U.table({ columns: ['Variant', 'Child', 'Father', 'Mother', 'Child DP / GQ', 'Allele balance', 'Quality'],
            rows: dn.map(function (r) { return [r.chrom + ':' + r.pos + ' ' + r.ref + '>' + r.alt, r.child_gt, r.father_gt,
              r.mother_gt, r.child_dp + ' / ' + r.child_gq, r.child_ab, r.quality]; }) }, 8));
      } else {
        add(d, el('p', { 'class': 'dim', 'data-trio': 'de_novo' }, empty(c.child, 'de_novo', 'de_novo.tsv')));
      }
      var ch = (f.compound_het || []).filter(function (r) { return r.child === c.child; });
      add(d, el('h3', {}, 'Compound-het pairs'));
      if (ch.length) {
        add(d, el('p', { 'class': 'dim' }, 'Two candidate variants in one gene, in trans first: one from each parent.'));
          add(d, U.table({ columns: ['Gene', 'Phase', 'Variant A', 'From', 'Variant B', 'From'],
            rows: ch.map(function (r) { return [r.gene, r.phase, r.variant_a, r.a_from, r.variant_b, r.b_from]; }) }, 8));
      } else {
        add(d, el('p', { 'class': 'dim', 'data-trio': 'compound_het' }, empty(c.child, 'compound_het', 'compound_het.tsv')));
      }
      add(s, d);
    });
    return s;
  }

  // The members side by side: what each one's report says of the genes that differ most between people.
  function sideBySide(f) {
    var sb = f.side_by_side;
    var s = section('Side by side', 'side-by-side');
    if (!sb || !sb.rows.length) { add(s, el('p', { 'class': 'dim' }, 'No member has a summary yet, or no summary lists these results.')); return s; }
    add(s, el('p', { 'class': 'dim' }, 'Each member’s own report, one column per person: CYP2D6, the CPIC genes ' +
      '(phenotypes in CPIC’s words), HLA, the haplogroups and the blood groups.'));
    // A value from a section older than its input (state "stale") says so in its own cell.
    var group = null, rows = [], old = 0;
    sb.rows.forEach(function (r) {
      rows.push([r.group === group ? '' : r.group, r.label].concat(r.values.map(function (v, i) {
        var m = sb.members[i];
        if (!(r.stale || [])[i]) { return v; }
        old++;
        return v + STALE;
      })));
      group = r.group;
    });
    add(s, U.table({ columns: ['', 'Gene or system'].concat(sb.members), rows: rows }, 14));
    if (old) {
      add(s, el('p', { 'class': 'dim', 'data-stale': String(old) }, 'Marked' + STALE + ': the member’s result file is from an ' +
        'earlier run, and the latest run did not write it again; that member’s page says why.'));
    }
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
    // A redraw (a result opened, a child's results answered) keeps each child's fold as the reader left it.
    var folds = {};
    Array.prototype.forEach.call(main.querySelectorAll('details[data-child]'), function (x) { folds[x.getAttribute('data-child')] = x.open; });
    var y = window.scrollY;
    main.textContent = '';
    var page = el('section', { 'class': 'family-page', 'data-plugin': 'family', 'aria-label': 'Family ' + f.state.name });
    if (f.state.error) {
      add(page, el('p', {}, f.state.error));
      add(main, page);
      root.setAttribute('data-pgp-status', 'ready');
      return;
    }
    var cover = el('div', { 'class': 'famcover' });
    var cl = el('div', {});
    add(cl, el('div', { 'class': 'kicker' }, 'Family report \u00b7 in five chapters'));
    add(cl, el('h1', { 'class': 'h-display' }, 'Family ' + f.state.name));
    var lp = el('p', { 'class': 'lede dropcap', id: 'family-lead' });
    lead(f).forEach(function (x) { add(lp, x[1] ? el('b', {}, x[0]) : x[0]); });
    add(cl, lp);
    var toc = el('ol', { 'class': 'famtoc' });
    [['ch-ped', 'I', 'Pedigree and members'], ['ch-related', 'II', 'Who is related how'], ['ch-couples', 'III', 'Couple carrier screen'],
      ['ch-children', 'IV', 'Each child: what they inherited'], ['ch-side-by-side', 'V', 'Side by side']].forEach(function (x) {
      add(toc, add(el('li', {}), add(el('a', { href: '#' + x[0] }), el('span', { 'class': 'secnum' }, x[1]), el('span', {}, x[2]))));
    });
    add(cover, cl, add(el('aside', {}), el('div', { 'class': 'fig-label' }, 'The chapters'), toc));
    add(page, cover);
    var top = el('div', { 'class': 'family-top' });
    var ped = section('Pedigree', 'ped');
    add(ped, pedigree(f), analyses(f));
    add(top, ped, members(f));
    add(page, top);
    add(page, related(f), couples(f), children(f), sideBySide(f));
    if ((f.state.plan || []).length || f.state.reports.length) { add(page, more(f)); }
    add(main, page);
    Array.prototype.forEach.call(main.querySelectorAll('details[data-child]'), function (x) {
      var k = x.getAttribute('data-child');
      if (k in folds) { x.open = folds[k]; }
    });
    window.scrollTo(0, y);
    // What the page knows of each child's results (for the tests): "DEMO-003:known".
    root.setAttribute('data-pgp-kids', U.familyStatus());
    root.setAttribute('data-pgp-status', 'ready');
  }


  // The page asks again every 30 s and redraws only when family.json changed (a member's run moved on,
  // step 39 finished). A failed first request says so, with a retry; a later one keeps what is shown.
  var drawnKey = null;
  function failed() {
    if (drawnKey !== null) { return; }
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
      // Each child's results are asked again with the family: a result to act on can appear while a run goes.
      if (k !== drawnKey) { drawnKey = k; draw(f); }
    }).catch(failed);
  }
  function start() {
    fetchFamily();
    window.setInterval(function () { if (!document.hidden) { fetchFamily(); } }, 30000);
  }
  if (document.readyState === 'loading') { document.addEventListener('DOMContentLoaded', start); } else { start(); }
  });
}());
