/* KINFOLK shared helpers: URL state, kinase index, search, mobile menu. */
(function () {
  "use strict";

  const K = (window.KINFOLK = window.KINFOLK || {});

  // ---------------------------------------------------------------- URL state
  K.state = {
    get(key, fallback) {
      const v = new URLSearchParams(location.search).get(key);
      return v === null ? fallback : v;
    },
    set(changes, replace = true) {
      const p = new URLSearchParams(location.search);
      for (const [k, v] of Object.entries(changes)) {
        if (v === null || v === undefined || v === "") p.delete(k);
        else p.set(k, v);
      }
      const url = `${location.pathname}${p.toString() ? "?" + p : ""}${location.hash}`;
      (replace ? history.replaceState : history.pushState).call(history, null, "", url);
    },
  };

  // ---------------------------------------------------------------- data
  let kinasePromise = null;
  K.kinases = function () {
    if (!kinasePromise) {
      kinasePromise = fetch("/api/kinases").then((r) => r.json()).then((rows) => {
        K.byId = new Map(rows.map((r) => [r.domain_id, r]));
        return rows;
      });
    }
    return kinasePromise;
  };
  K.json = (url) => fetch(url).then((r) => (r.ok ? r.json() : Promise.reject(new Error(`${r.status} ${url}`))));

  K.escape = (s) =>
    String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);

  K.fmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? "n/a" : Number(v).toFixed(d));

  // Marker glyphs (shape carries meaning as well as colour).
  K.colours = {
    agree: "#3F6B52", partial: "#C98F2B", disagree: "#A9482A", none: "none", pending: "#C9D0C1",
    esm: "#E07A55", afdb: "#C98F2B", pdb: "#1F3326", ink: "#1F3326", rule: "#C9D0C1",
  };
  K.agreementLabel = { agree: "agrees", partial: "partial", disagree: "disagrees", none: "no PDB structure", pending: "pending", missing: "no model" };

  // ---------------------------------------------------------------- search
  function initSearch() {
    const input = document.getElementById("kinase-search");
    const list = document.getElementById("search-results");
    if (!input || !list) return;
    let items = [];
    let active = -1;

    function render(q) {
      const needle = q.trim().toLowerCase();
      list.innerHTML = "";
      active = -1;
      if (!needle) {
        input.setAttribute("aria-expanded", "false");
        return;
      }
      K.kinases().then((rows) => {
        items = rows
          .filter((r) =>
            [r.domain_id, r.gene, r.uniprot, r.family, r.grp].some((f) => String(f || "").toLowerCase().includes(needle)))
          .slice(0, 12);
        list.innerHTML = items
          .map((r, i) => `<li><a id="sr-${i}" role="option" href="/k/${encodeURIComponent(r.domain_id)}" aria-selected="false">
              <span>${K.escape(r.domain_id)}</span><span class="mono" style="opacity:.8">${K.escape(r.grp)}</span></a></li>`)
          .join("");
        input.setAttribute("aria-expanded", items.length ? "true" : "false");
      });
    }

    function move(delta) {
      const opts = list.querySelectorAll("a");
      if (!opts.length) return;
      if (active >= 0) opts[active].setAttribute("aria-selected", "false");
      active = (active + delta + opts.length) % opts.length;
      opts[active].setAttribute("aria-selected", "true");
      input.setAttribute("aria-activedescendant", opts[active].id);
    }

    input.addEventListener("input", () => render(input.value));
    input.addEventListener("keydown", (e) => {
      if (e.key === "ArrowDown") { e.preventDefault(); move(1); }
      else if (e.key === "ArrowUp") { e.preventDefault(); move(-1); }
      else if (e.key === "Enter") {
        const opts = list.querySelectorAll("a");
        const target = active >= 0 ? opts[active] : opts[0];
        if (target) { e.preventDefault(); location.href = target.href; }
      } else if (e.key === "Escape") { input.value = ""; render(""); }
    });
  }

  function initMenu() {
    const btn = document.getElementById("menu-toggle");
    const bar = document.getElementById("sidebar");
    if (!btn || !bar) return;
    btn.addEventListener("click", () => {
      const open = bar.dataset.open !== "true";
      bar.dataset.open = String(open);
      btn.setAttribute("aria-expanded", String(open));
    });
  }

  // Mol* inserts an empty logo link; give it an accessible name wherever a viewer appears.
  function labelMolstar(root) {
    root.querySelectorAll("a.msp-logo:not([aria-label])").forEach((a) => {
      a.setAttribute("aria-label", "Mol* viewer website (opens in a new tab)");
      a.setAttribute("rel", "noopener noreferrer");
    });
  }

  document.addEventListener("DOMContentLoaded", () => {
    initSearch();
    initMenu();
    labelMolstar(document);
    new MutationObserver(() => labelMolstar(document)).observe(document.body, { childList: true, subtree: true });
  });
})();
