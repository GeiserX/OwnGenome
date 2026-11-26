// Browser-side property helper. Walk unknown fields too; callers mark only immutable controls/ownership.
(function () {
  'use strict';
  function walk(value, keep, path) {
    path = path || [];
    if (keep(path, value)) { return JSON.parse(JSON.stringify(value)); }
    if (Array.isArray(value)) {
      var items = value.map(function (v, i) { return walk(v, keep, path.concat(i)); });
      if (items.length) { items.push(walk(value[0], keep, path.concat(value.length))); }
      else { items.push('future-array-value'); }
      return items;
    }
    if (value && typeof value === 'object') {
      var out = { future_payload_field: 'future-value' };
      Object.keys(value).reverse().forEach(function (k) { out[k] = walk(value[k], keep, path.concat(k)); });
      return out;
    }
    if (typeof value === 'string') { return value + '_CHANGED'; }
    if (typeof value === 'number') { return value + 37; }
    if (typeof value === 'boolean') { return !value; }
    return 'changed-null';
  }
  function snapshot() {
    // Attribute order is part of this snapshot. No fields are removed from the gated document.
    return {
      dom: document.documentElement.outerHTML,
      text: document.body.innerText,
      title: document.title,
      address: location.pathname + location.search + location.hash,
      canvas: Array.from(document.querySelectorAll('canvas')).map(function (c) { return c.toDataURL(); }),
      svg: Array.from(document.querySelectorAll('svg')).map(function (s) { return [s.outerHTML, s.clientWidth, s.clientHeight]; }),
      storage: [Object.assign({}, localStorage), Object.assign({}, sessionStorage)],
      requests: performance.getEntriesByType('resource').map(function (r) { return new URL(r.name).pathname; }).sort()
    };
  }
  window.GateProperty = { walk: walk, snapshot: snapshot };
}());
