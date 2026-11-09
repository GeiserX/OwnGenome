// live.js: the run's progress on a person's page (assets/live/index.html, served by bin/owngenome_live.py):
// the step timeline (one line per step), the step detail, the wire of what each step found, the event
// stream and the replay, in the report's "How the run went" section. One data source feeds the page: in live mode state.json plus the run's event stream;
// with <html data-mode="replay"> a frames.json of recorded states played on a timer. person.js draws
// the results from the same states. No outside request, no inline script, no eval. Plugins call
// window.PGP.register({name, onState, onEvent, panel, slot}); a panel goes into the element whose id
// is slot when the page has one (the search box in the header), else beside the step timeline.
(function () {
  'use strict';
  var root = document.documentElement;
  var blocked = [];
  document.addEventListener('securitypolicyviolation', function (e) {
    blocked.push(e.violatedDirective + ' ' + (e.blockedURI || ''));
    root.setAttribute('data-pgp-blocked', blocked.join(' | '));
  });

  // One table maps each run_state state to its class and the words the legend shows.
  var STATES = {
    'waiting': ['s-waiting', 'waiting'],
    'running': ['s-running', 'running'],
    'done': ['s-done', 'done'],
    'failed': ['s-failed', 'failed'],
    'skipped': ['s-skipped', 'skipped'],
    'not finished': ['s-notfinished', 'not finished'],
    'not planned': ['s-notplanned', 'not planned']
  };
  var SVGNS = 'http://www.w3.org/2000/svg';
  var plugins = [];
  var last = null, selected = null, boxes = {};

  function el(tag, attrs, text) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    if (text !== undefined && text !== null) { e.textContent = String(text); }
    return e;
  }
  function svg(tag, attrs) {
    var e = document.createElementNS(SVGNS, tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    return e;
  }
  function stateOf(s) { return STATES[s] || ['s-other', s]; }

  function shown(steps) {
    var v = steps.filter(function (s) { return s.state !== 'not planned'; });
    return v.length ? v : steps;
  }

  // The run as a timeline, one line per step (the Journal's step timeline): a step's bar runs from its
  // start to its end; a running step's bar ends at now, with its usual hours drawn as an outline after
  // it; a step still waiting has no bar yet. Steps are ordered by start, the waiting ones last, in the
  // order assets/steps.tsv plans them. Each line is a button (g.box) that opens the step's detail.
  function ms(t) { var v = t ? Date.parse(t) : NaN; return isFinite(v) ? v : null; }
  function clock(v) { var d = new Date(v); return ('0' + d.getUTCHours()).slice(-2) + ':' + ('0' + d.getUTCMinutes()).slice(-2); }
  function timeline(steps, nowText, running) {
    // A finished run ends at its last step, not at the moment the page was asked for.
    var now = running ? ms(nowText) : null;
    var known = steps.map(function (s) { return ms(s.started_utc); }).filter(function (v) { return v !== null; });
    var t0 = known.length ? Math.min.apply(null, known) : null;
    var ends = steps.map(function (s) { return ms(s.ended_utc) || (s.state === 'running' ? now : null); })
      .filter(function (v) { return v !== null; });
    var t1 = Math.max.apply(null, ends.concat(known).concat(now !== null ? [now] : []).concat(t0 !== null ? [t0 + 3600e3] : [0]));
    var order = steps.map(function (s, i) { return { s: s, i: i, t: ms(s.started_utc) }; });
    order.sort(function (a, b) {
      if (a.t === null && b.t === null) { return a.i - b.i; }
      if (a.t === null) { return 1; } if (b.t === null) { return -1; }
      return a.t - b.t || a.i - b.i;
    });
    // Leave room after now for the usual hours of what still runs.
    var tail = 0;
    steps.forEach(function (s) {
      if (s.state === 'running' && s.hours && ms(s.started_utc) !== null) { tail = Math.max(tail, ms(s.started_utc) + s.hours * 3600e3 - t1); }
    });
    return { order: order.map(function (o) { return o.s; }), t0: t0, t1: t1 + Math.max(0, Math.min(tail, (t1 - (t0 || t1)) * 0.6)), now: now };
  }
  var LABEL_W = 190, ROW_H = 16, TOP = 8, AXIS_H = 26;

  function drawMap(steps) {
    var map = document.getElementById('map');
    var T = timeline(steps, last && (last.generated_utc || null), !!(last && last.run && last.run.state === 'running'));
    var h = TOP + T.order.length * ROW_H + AXIS_H;
    // The chart takes the width the page gives it (at least 300 px of time), so no step sits off screen.
    var CHART_W = Math.max(300, (map.clientWidth || 962) - LABEL_W - 12);
    var W = LABEL_W + CHART_W + 12;
    var s = svg('svg', { viewBox: '0 0 ' + W + ' ' + h, width: W, height: h, role: 'group', 'aria-label': 'Steps of the run' });
    var span = T.t0 === null ? 1 : Math.max(T.t1 - T.t0, 60e3);
    function x(v) { return LABEL_W + 8 + (T.t0 === null ? 0 : (v - T.t0) / span * (CHART_W - 8)); }
    var axis = svg('g', { 'class': 'axis' });
    if (T.t0 !== null) {
      var stepH = span > 36 * 3600e3 ? 6 : span > 12 * 3600e3 ? 3 : span > 4 * 3600e3 ? 1 : 0.5;
      for (var t = Math.ceil(T.t0 / (stepH * 3600e3)) * stepH * 3600e3; t <= T.t1; t += stepH * 3600e3) {
        axis.appendChild(svg('line', { x1: x(t), x2: x(t), y1: TOP - 4, y2: h - AXIS_H + 2 }));
        // No hour label under the "now" label of a running run.
        if (T.now !== null && last && last.run && last.run.state === 'running' && Math.abs(x(t) - x(T.now)) < 56) { continue; }
        var lab = svg('text', { x: x(t), y: h - 8, 'text-anchor': 'middle' });
        lab.textContent = clock(t);
        axis.appendChild(lab);
      }
    }
    s.appendChild(axis);
    boxes = {};
    T.order.forEach(function (st, i) {
      var y = TOP + i * ROW_H;
      var g = svg('g', { 'class': 'box', 'data-step': st.step, tabindex: '0', role: 'button' });
      g.appendChild(svg('rect', { 'class': 'hit', x: 0, y: y, width: W, height: ROW_H }));
      var t = svg('text', { 'class': 'lbl', x: LABEL_W, y: y + ROW_H - 4, 'text-anchor': 'end' });
      var label = st.label || st.step;
      t.textContent = label.length > 30 ? label.slice(0, 29) + '…' : label;
      g.appendChild(t);
      var a = ms(st.started_utc), b = ms(st.ended_utc) || (st.state === 'running' ? T.now : null);
      if (a !== null && T.t0 !== null) {
        var x1 = x(a), x2 = Math.max(x(b !== null ? b : a), x1 + 3);
        g.appendChild(svg('rect', { 'class': 'bar', x: x1, y: y + 4, width: x2 - x1, height: ROW_H - 8, rx: 1.5 }));
        if (st.state === 'running' && st.hours) {
          var u = x(a + st.hours * 3600e3);
          if (u > x2) { g.appendChild(svg('rect', { 'class': 'usual', x: x2, y: y + 4.5, width: Math.min(u, W - 4) - x2, height: ROW_H - 9, rx: 1.5 })); }
        }
      } else {
        g.appendChild(svg('rect', { 'class': 'bar', x: LABEL_W + 8, y: y + 5, width: 0, height: 0 }));
      }
      g.appendChild(svg('title', {}));
      g.addEventListener('click', function () { selected = st.step; detail(); });
      g.addEventListener('keydown', function (e) { if (e.key === 'Enter' || e.key === ' ') { selected = st.step; detail(); } });
      s.appendChild(g);
      boxes[st.step] = g;
    });
    if (T.now !== null && T.t0 !== null && last && last.run && last.run.state === 'running') {
      s.appendChild(svg('line', { 'class': 'nowline', x1: x(T.now), x2: x(T.now), y1: TOP - 4, y2: h - AXIS_H + 2 }));
      var nl = svg('text', { 'class': 'nowlab', x: x(T.now), y: h - 8, 'text-anchor': 'middle' });
      nl.textContent = 'now ' + clock(T.now);
      s.appendChild(nl);
    }
    map.textContent = '';
    map.appendChild(s);
    steps.forEach(function (st) { paint(st.step, st.state, st.label); });
  }

  function paint(step, state, label) {
    var g = boxes[step];
    if (!g) { return; }
    var cls = stateOf(state);
    g.setAttribute('class', 'box ' + cls[0] + (step === selected ? ' selected' : ''));
    g.setAttribute('data-state', state);
    var title = g.querySelector('title');
    if (title && label) { title.textContent = step + ' ' + label + ': ' + cls[1]; }
    if (label) { g.setAttribute('aria-label', step + ' ' + label + ': ' + cls[1]); }
  }

  function legend() {
    var l = document.getElementById('legend');
    if (l.childNodes.length) { return; }
    Object.keys(STATES).forEach(function (k) {
      var i = el('span', { 'class': 'key ' + STATES[k][0] });
      i.appendChild(el('span', { 'class': 'swatch' }));
      i.appendChild(document.createTextNode(STATES[k][1]));
      l.appendChild(i);
    });
  }

  function detail() {
    var box = document.getElementById('detail');
    if (!last || !selected) { return; }
    var st = last.steps.filter(function (s) { return s.step === selected; })[0];
    if (!st) { return; }
    box.textContent = '';
    box.appendChild(el('h2', {}, 'Step ' + st.step + ': ' + st.label));
    var dl = el('dl', {});
    function row(k, v) { if (v === null || v === undefined || v === '') { return; } dl.appendChild(el('dt', {}, k)); dl.appendChild(el('dd', {}, v)); }
    row('State', stateOf(st.state)[1] + (st.reason ? ' (' + st.reason + ')' : ''));
    row('Note', st.note);
    row('Started', st.started_utc);
    row('Ended', st.ended_utc);
    row('Usual hours', st.hours);
    box.appendChild(dl);
    var procs = Object.keys(st.processes || {});
    if (procs.length) {
      var t = el('table', { 'class': 'tasks' });
      var h = el('tr', {});
      ['Process', 'waiting', 'running', 'done', 'cached', 'failed'].forEach(function (c) { h.appendChild(el('th', {}, c)); });
      t.appendChild(h);
      procs.forEach(function (p) {
        var r = el('tr', {});
        r.appendChild(el('td', {}, p));
        ['waiting', 'running', 'done', 'cached', 'failed'].forEach(function (c) { r.appendChild(el('td', {}, st.processes[p][c])); });
        t.appendChild(r);
      });
      box.appendChild(t);
    }
    (st.sections || []).forEach(function (k) {
      var hl = (last.headlines || {})[k];
      if (hl) { box.appendChild(el('p', { 'class': 'headline tone-' + hl.tone }, hl.text)); }
      if (root.getAttribute('data-mode') === 'live') {
        var anchor = (last.anchors || {})[k] || k;
        box.appendChild(el('a', { href: 'report.html#' + encodeURIComponent(anchor) }, 'Open in the report: ' + k));
      }
    });
    // Repaint, never redraw: the focused step keeps its focus and the plugins' marks stay on.
    shown(last.steps).forEach(function (s) { paint(s.step, s.state, s.label); });
  }

  // The wire: what the steps said, newest first. Only an event that says something gets its own line (a
  // finding written, a warning, a failure); the plain "finished, see the report" ones fold into one line
  // that names the newest. Each line is cut to its first sentence; the step's detail has the whole text.
  // The newest WIRE_SHOWN lines show, the rest sit under "Every event".
  var WIRE_SHOWN = 8, WIRE_CHARS = 140;
  var PLAIN = /^(.*?):? finished, see the report\.?$/;
  function sentence(t) {
    var m = /^(.+?[.!?])(\s+[A-Z(]|$)/.exec(t || '');
    var s = m ? m[1] : (t || '');
    return s.length > WIRE_CHARS ? s.slice(0, WIRE_CHARS - 1).replace(/\s+\S*$/, '') + '…' : s;
  }
  function wireLine(ev) {
    var li = el('li', { 'class': 'tone-' + ev.tone, 'data-section': ev.key || '', tabindex: '0', role: 'button', title: ev.full });
    li.appendChild(el('span', { 'class': 'step' }, ev.t ? clock(ev.t) : ev.step));
    li.appendChild(document.createTextNode(' ' + ev.text));
    function open() { selected = ev.step; detail(); }
    li.addEventListener('click', open);
    li.addEventListener('keydown', function (e) { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(); } });
    return li;
  }
  function found(state) {
    var ul = document.getElementById('found');
    var hl = state.headlines || {};
    var said = [], plain = [];
    state.steps.forEach(function (s, i) {
      var t = ms(s.ended_utc);
      (s.sections || []).forEach(function (k) {
        if (!hl[k]) { return; }
        var m = PLAIN.exec(hl[k].text || '');
        var ev = { key: k, step: hl[k].step || s.step, tone: hl[k].tone, t: t, i: i, full: hl[k].text };
        if (m && hl[k].tone !== 'attention') { ev.name = m[1] || s.label; plain.push(ev); }
        else { ev.text = sentence(hl[k].text); said.push(ev); }
      });
      if (s.state === 'failed') { said.push({ step: s.step, tone: 'attention', t: t, i: i, text: (s.label || s.step) + ' failed', full: s.reason || '' }); }
    });
    function newest(a, b) { return (b.t || 0) - (a.t || 0) || b.i - a.i; }
    plain.sort(newest);
    if (plain.length) {
      var first = plain[0], since = plain.map(function (e) { return e.t; }).filter(Boolean).sort()[0];
      said.push({ key: first.key, step: first.step, tone: 'info', t: first.t, i: first.i, full: first.full,
        text: first.name + (plain.length > 1 ? ' and ' + (plain.length - 1) + ' more step' + (plain.length > 2 ? 's' : '') +
          ' finished' + (since ? ' since ' + clock(since) : '') : ' finished') + ', see the report' });
    }
    said.sort(newest);
    ul.textContent = '';
    if (!said.length) { ul.appendChild(el('li', { 'class': 'muted' }, 'Nothing finished yet.')); return; }
    said.slice(0, WIRE_SHOWN).forEach(function (ev) { ul.appendChild(wireLine(ev)); });
    var rest = said.slice(WIRE_SHOWN);
    if (rest.length || plain.length > 1) {
      var li = el('li', { 'class': 'wire-more' });
      var d = el('details', {});
      d.appendChild(el('summary', {}, 'Every event'));
      var all = el('ul', { 'class': 'found' });
      rest.forEach(function (ev) { all.appendChild(wireLine(ev)); });
      plain.forEach(function (ev) { all.appendChild(wireLine({ key: ev.key, step: ev.step, tone: 'info', t: ev.t, text: ev.name + ' finished', full: ev.full })); });
      d.appendChild(all);
      li.appendChild(d);
      ul.appendChild(li);
    }
  }

  function show(state) {
    last = state;
    window.PGP.state = state;
    legend();
    var badge = document.getElementById('run-state');
    badge.textContent = state.run.state;
    badge.setAttribute('class', 'badge run-' + state.run.state.replace(/\s+/g, '-'));
    drawMap(shown(state.steps));
    found(state);
    detail();
    plugins.forEach(function (p) { if (p.onState) { try { p.onState(state); } catch (e) { console.error(p.name, e); } } });
    root.setAttribute('data-pgp-status', 'ready');
  }

  // --- live: state.json and the event stream ---
  var timer = null;
  function refresh() {
    fetch('state.json', { credentials: 'same-origin', cache: 'no-store' }).then(function (r) {
      if (r.status === 401) { window.location.href = '/login?next=' + encodeURIComponent(window.location.pathname); return null; }
      return r.ok ? r.json() : null;
    }).then(function (s) { if (s) { show(s); } }).catch(function () { root.setAttribute('data-pgp-status', 'offline'); });
  }
  function later() { if (timer) { clearTimeout(timer); } timer = setTimeout(refresh, 800); }
  var opens = 0;
  function live() {
    refresh();
    stream();
  }
  function stream() {
    var es = new EventSource('events');
    root.setAttribute('data-pgp-streams', String(++opens));
    es.onmessage = function (e) {
      var m;
      try { m = JSON.parse(e.data); } catch (x) { return; }
      // Paint at once what a new message says, then let state.json settle it. A message the
      // last state.json already counted (the stream replays from the start) is not painted.
      var fresh = !last || Number(e.lastEventId) > (last.events || 0);
      if (fresh && m.step && (m.event === 'process_started' || m.event === 'process_submitted' || m.state === 'started')) { paint(m.step, 'running'); }
      plugins.forEach(function (p) { if (p.onEvent) { try { p.onEvent(m); } catch (x) { console.error(p.name, x); } } });
      later();
    };
    // EventSource retries a dropped stream by itself, but any refusal (a 401 when the session ended,
    // a 503 when every stream slot is taken) closes it for good. Then ask state.json, whose 401 goes
    // to the login, and open the stream again a little later.
    es.onerror = function () {
      if (es.readyState === EventSource.CLOSED) {
        root.setAttribute('data-pgp-status', 'offline');
        refresh();
        setTimeout(stream, 3000);
      }
    };
  }

  // --- replay: frames.json on a timer ---
  function replay() {
    var bar = document.getElementById('replay');
    bar.hidden = false;
    var report = document.getElementById('report-link');
    if (report) { report.hidden = true; }
    fetch('frames.json').then(function (r) { return r.json(); }).then(function (data) {
      var frames = data.frames || [];
      var slider = document.getElementById('replay-slider'), play = document.getElementById('replay-play');
      var at = 0, tick = null;
      slider.max = String(Math.max(frames.length - 1, 0));
      function go(i) {
        at = Math.min(Math.max(i, 0), frames.length - 1);
        slider.value = String(at);
        document.getElementById('replay-time').textContent = frames[at].t || '';
        show(frames[at].state);
        root.setAttribute('data-pgp-frame', String(at));
      }
      function stop() { if (tick) { clearInterval(tick); tick = null; } play.textContent = 'Play'; }
      play.addEventListener('click', function () {
        if (tick) { stop(); return; }
        if (at >= frames.length - 1) { go(0); }
        play.textContent = 'Pause';
        tick = setInterval(function () { if (at >= frames.length - 1) { stop(); } else { go(at + 1); } }, Number(data.interval_ms) || 1000);
      });
      slider.addEventListener('input', function () { stop(); go(Number(slider.value)); });
      if (frames.length) { go(0); }
      if (data.autoplay) { play.click(); }
    });
  }

  window.PGP = {
    state: null,
    register: function (p) {
      plugins.push(p);
      if (p.panel) {
        var box = el('section', { 'class': 'card plugin', 'data-plugin': p.name || '' });
        var slot = (p.slot && document.getElementById(p.slot)) || document.getElementById('plugin-panels');
        slot.appendChild(box);
        try { p.panel(box); } catch (e) { console.error(p.name, e); }
      }
      if (last && p.onState) { p.onState(last); }
    }
  };

  // The timeline is drawn to the width the page gives it: a fold opened after the first draw draws it again.
  var fold = document.getElementById('progress-fold');
  if (fold) { fold.addEventListener('toggle', function () { if (fold.open && last) { drawMap(shown(last.steps)); detail(); } }); }

  function startPage() {
    if (root.getAttribute('data-mode') === 'replay') { replay(); } else { live(); }
  }
  if (document.readyState === 'loading') { document.addEventListener('DOMContentLoaded', startPage); } else { startPage(); }
}());
