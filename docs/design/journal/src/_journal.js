/* OwnGenome "Journal": tiny shared helpers (tooltip, reveals, inline figures). Synthetic data only. */
(function () {
  const NS = 'http://www.w3.org/2000/svg';
  const CHR = [['1',248.96,123.4],['2',242.19,93.9],['3',198.3,90.9],['4',190.21,50],['5',181.54,48.8],['6',170.81,59.8],['7',159.35,60.1],['8',145.14,45.2],['9',138.39,43],['10',133.8,39.8],['11',135.09,53.4],['12',133.28,35.5],['13',114.36,17.7],['14',107.04,17.2],['15',101.99,19],['16',90.34,36.8],['17',83.26,25.1],['18',80.37,18.5],['19',58.62,26.2],['20',64.44,28.1],['21',46.71,12],['22',50.82,15],['X',156.04,60.6],['Y',57.23,10.4]];
  const LV = { act: '#d03b3b', med: '#ec835a', rev: '#fab219', car: '#9085e9', info: '#c3c2b7' };
  window.OG = { CHR, LV };

  function el(tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  OG.el = el;

  // deterministic PRNG for synthetic paintings
  OG.rng = function (seed) { let s = seed >>> 0; return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296); };

  /* Karyogram: vertical chromosomes in one row, findings as level glyphs beside them.
     marks: [{chr:'13', mb:32.3, lv:'act', label:'BRCA2', tip:'...'}] */
  OG.karyogram = function (host, marks, o = {}) {
    const chrs = CHR.filter(c => !(o.noY && c[0] === 'Y'));
    const W = o.width || 640, H = o.height || 170, top = 22, bottom = 22;
    const slot = (W - (o.padRight ?? 56)) / chrs.length, maxL = 249, scale = (H - top - bottom) / maxL;
    const s = el('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', role: 'img', 'aria-label': o.aria || 'Chromosome map of findings' }, host);
    chrs.forEach((c, i) => {
      const x = i * slot + slot / 2 - 3.5, h = c[1] * scale, cy = top + c[2] * scale;
      el('rect', { x, y: top, width: 7, height: c[2] * scale - 1.2, rx: 3.5, fill: '#2c2d33' }, s);
      el('rect', { x, y: cy + 1.2, width: 7, height: h - c[2] * scale - 1.2, rx: 3.5, fill: '#2c2d33' }, s);
      const t = el('text', { x: x + 3.5, y: H - 6, 'text-anchor': 'middle', 'font-size': 10, fill: '#8f8d86' }, s); t.textContent = c[0];
    });
    marks.forEach(m => {
      const i = chrs.findIndex(c => c[0] === m.chr); if (i < 0) return;
      const x = i * slot + slot / 2 + 7, y = top + m.mb * scale;
      const g = el('g', { 'data-tip': m.tip || m.label }, s);
      el('line', { x1: x - 10.5, x2: x - 1, y1: y, y2: y, stroke: LV[m.lv], 'stroke-width': 2 }, g);
      glyph(g, m.lv, x + 4, y);
      if (m.label && !o.noLabels) { const t = el('text', { x: x + 11, y: y + 3.5, 'font-size': 11, 'font-weight': 600, fill: '#f0efec', stroke: o.bg || '#0e0f10', 'stroke-width': 4, 'paint-order': 'stroke', 'stroke-linejoin': 'round' }, g); t.textContent = m.label; }
      el('rect', { x: x - 12, y: y - 9, width: 60, height: 18, fill: 'transparent' }, g);
    });
    return s;
  };
  function glyph(g, lv, x, y) {
    const c = LV[lv];
    if (lv === 'act') el('rect', { x: x - 4, y: y - 4, width: 8, height: 8, fill: c, transform: `rotate(45 ${x} ${y})` }, g);
    else if (lv === 'med') el('rect', { x: x - 5.5, y: y - 3, width: 11, height: 6, rx: 3, fill: c }, g);
    else if (lv === 'rev') el('path', { d: `M${x} ${y - 5}L${x + 5} ${y + 4}L${x - 5} ${y + 4}Z`, fill: c }, g);
    else if (lv === 'car') { el('circle', { cx: x - 2.6, cy: y, r: 3.4, fill: 'none', stroke: c, 'stroke-width': 1.6 }, g); el('circle', { cx: x + 2.6, cy: y, r: 3.4, fill: 'none', stroke: c, 'stroke-width': 1.6 }, g); }
    else el('circle', { cx: x, cy: y, r: 4, fill: 'none', stroke: c, 'stroke-width': 1.5 }, g);
  }
  OG.glyph = glyph;

  /* Percentile strip: 0-100, dot at value. No zones (a percentile is not a risk class). */
  OG.pstrip = function (host, value, o = {}) {
    const W = o.width || 260, H = 22;
    const s = el('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H }, host);
    el('line', { x1: 0, x2: W, y1: 11, y2: 11, stroke: 'rgba(255,255,255,.16)', 'stroke-width': 1 }, s);
    [0, 25, 50, 75, 100].forEach(t => el('line', { x1: t / 100 * (W - 2) + 1, x2: t / 100 * (W - 2) + 1, y1: 7, y2: 15, stroke: 'rgba(255,255,255,.16)' }, s));
    const x = value / 100 * (W - 12) + 6;
    el('circle', { cx: x, cy: 11, r: 5.5, fill: o.color || '#f0efec', stroke: '#15161a', 'stroke-width': 2, 'data-tip': o.tip || `${value}th percentile` }, s);
    return s;
  };

  /* Sparkline */
  OG.spark = function (host, data, o = {}) {
    const W = o.width || 120, H = o.height || 28, max = o.max || Math.max(...data), min = o.min ?? 0;
    const s = el('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H }, host);
    const pts = data.map((v, i) => [i / (data.length - 1) * (W - 4) + 2, H - 3 - (v - min) / (max - min) * (H - 6)]);
    el('path', { d: 'M' + pts.map(p => p.join(' ')).join('L') + `L${W - 2} ${H}L2 ${H}Z`, fill: o.fill || 'rgba(240,239,236,.06)' }, s);
    el('path', { d: 'M' + pts.map(p => p.join(' ')).join('L'), fill: 'none', stroke: o.color || '#c3c2b7', 'stroke-width': 1.5, 'stroke-linejoin': 'round' }, s);
    const l = pts[pts.length - 1]; el('circle', { cx: l[0], cy: l[1], r: 2.6, fill: o.color || '#c3c2b7' }, s);
    return s;
  };

  /* Chromosome painting: 22 autosomes x 2 haplotypes, segments coloured by reference group. */
  OG.painting = function (host, groups, seed, o = {}) {
    const r = OG.rng(seed), W = o.width || 640, cols = o.cols || 1, rowH = 5, gap = 2, chrGap = 9, left = 22, colGap = 28;
    const auto = CHR.slice(0, 22), perCol = Math.ceil(auto.length / cols), colW = (W - colGap * (cols - 1)) / cols;
    const maxL = 249, sc = (colW - left - 4) / maxL;
    const H = perCol * (rowH * 2 + gap + chrGap);
    const s = el('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', role: 'img', 'aria-label': 'Chromosome painting' }, host);
    const cum = []; let acc = 0; groups.forEach(g => { acc += g.f; cum.push(acc); });
    const pick = () => { const v = r() * acc; return groups[cum.findIndex(c => v <= c)]; };
    auto.forEach((c, i) => {
      const col = Math.floor(i / perCol), ox = col * (colW + colGap);
      const y0 = (i % perCol) * (rowH * 2 + gap + chrGap) + 2;
      const t = el('text', { x: ox, y: y0 + rowH + 4, 'font-size': 10, fill: '#8f8d86' }, s); t.textContent = c[0];
      for (let hap = 0; hap < 2; hap++) {
        let pos = 0; const y = y0 + hap * (rowH + gap); const segs = [];
        while (pos < c[1]) {
          const g = pick(); const len = Math.min(c[1] - pos, 3 + r() * 22);
          const last = segs[segs.length - 1];
          if (last && last.g === g) last.len += len; else segs.push({ g, start: pos, len });
          pos += len;
        }
        segs.forEach((sg, k) => el('rect', { x: ox + left + sg.start * sc, y, width: Math.max(sg.len * sc - (k < segs.length - 1 ? 1.5 : 0), .8), height: rowH, rx: 1, fill: sg.g.c, 'data-tip': `chr${c[0]} copy ${hap + 1}: ${sg.g.n}` }, s));
      }
    });
    return s;
  };

  // tooltip + reveals
  document.addEventListener('DOMContentLoaded', () => {
    const tip = document.createElement('div'); tip.id = 'tip'; document.body.appendChild(tip);
    document.addEventListener('mousemove', e => {
      const t = e.target.closest && e.target.closest('[data-tip]');
      if (!t) { tip.style.opacity = 0; return; }
      tip.textContent = t.getAttribute('data-tip'); tip.style.opacity = 1;
      tip.style.left = Math.min(e.clientX + 14, innerWidth - 280) + 'px'; tip.style.top = (e.clientY + 16) + 'px';
    });
    document.querySelectorAll('[data-reveal]').forEach(b => b.addEventListener('click', () => {
      const tgt = document.getElementById(b.getAttribute('data-reveal'));
      if (tgt) { tgt.hidden = !tgt.hidden; if (b.dataset.once) b.closest('.veil')?.remove(); }
    }));
  });
})();
