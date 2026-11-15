// The page of scripts/owngenome-view.sh: reads config.json from the local server
// (assets/viewer/serve.py) and opens igv.js on it at the locus in the URL
// after '#', else the one the config names. The report's links set that
// part, so a click on a finding opens the reads at it.
//
// The state goes on <html data-pgp-status>: "loading", "ready" or
// "error: <why>"; a request the Content-Security-Policy blocked is added to
// data-pgp-blocked. tests/test_viewer.py reads both from a headless browser.
"use strict";

(function () {
  const root = document.documentElement;
  const blocked = [];

  document.addEventListener("securitypolicyviolation", function (e) {
    blocked.push(e.violatedDirective + " " + (e.blockedURI || ""));
    root.dataset.pgpBlocked = blocked.join(" | ");
  });

  function fail(why) {
    root.dataset.pgpStatus = "error: " + why;
    const box = document.getElementById("error");
    box.textContent = "The viewer could not start: " + why;
    box.style.display = "block";
  }

  function locusFromHash() {
    try {
      return decodeURIComponent(window.location.hash.slice(1)).trim();
    } catch (e) {
      return "";
    }
  }

  async function start() {
    if (typeof igv === "undefined") {
      fail("igv.min.js did not load");
      return;
    }
    const resp = await fetch("config.json", { cache: "no-store" });
    if (!resp.ok) {
      fail("config.json answered " + resp.status);
      return;
    }
    const cfg = await resp.json();
    cfg.tracks.forEach(function (track) {
      if (track.type === "alignment") {
        track.visibilityWindow = 50000;
      }
    });
    document.getElementById("sample").textContent = cfg.sample;
    document.title = cfg.sample + " - Genome viewer";
    const browser = await igv.createBrowser(document.getElementById("igv"), {
      reference: cfg.reference,
      locus: locusFromHash() || cfg.locus,
      tracks: cfg.tracks,
      showCenterGuide: true,
      // igv.js otherwise fetches its genome list from igv.org at start.
      loadDefaultGenomes: false,
    });
    function followLocus() {
      const locus = locusFromHash();
      if (locus) {
        browser.search(locus);
      }
    }
    // igv.js keeps bases per pixel on resize. While igv still shows the locus the URL asks
    // for, a resize applies that locus again. Once the reader pans, zooms or types a locus in
    // igv, igv keeps the reader's place, until the URL asks for a new locus. The events
    // come out of igv's shadow root, because pointer, wheel and key events are composed.
    let follow = true;
    ["pointerdown", "wheel", "keydown"].forEach(function (type) {
      document.getElementById("igv").addEventListener(type, function () {
        follow = false;
      }, { capture: true, passive: true });
    });
    window.addEventListener("hashchange", function () {
      follow = true;
      followLocus();
    });
    let resizeTimer;
    window.addEventListener("resize", function () {
      window.clearTimeout(resizeTimer);
      resizeTimer = window.setTimeout(function () {
        if (follow) {
          followLocus();
        }
      }, 100);
    });
    root.dataset.pgpTracks = String(cfg.tracks.length);
    root.dataset.pgpStatus = "ready";
  }

  start().catch(function (e) {
    fail(String((e && e.message) || e));
  });
})();
