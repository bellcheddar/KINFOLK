/* Plate V: pseudokinase gallery with flip cards and a shared Mol* dialog. */
(function () {
  "use strict";
  const K = window.KINFOLK;
  const CALL_CLASS = { restored: "disagree", partial: "partial", faithful: "agree" };
  let rows = [], viewer = null;

  function card(r) {
    const detail = JSON.parse(r.detail_json || "[]");
    const byMotif = {};
    detail.forEach((d) => ((byMotif[d.motif] = byMotif[d.motif] || []).push(d)));
    const calls = [["afdb", "■ AlphaFold DB", r.afdb_call], ["esm2", "▲ ESMFold2", r.esm_call], ["pdb", "○ PDB", r.pdb_call]]
      .map(([, label, c]) => `<li>${label}: ${c ? `<span class="chip ${CALL_CLASS[c] || ""}">${K.escape(c)}</span>` : '<span class="muted">n/a</span>'}</li>`).join("");
    const errs = Object.entries(byMotif).map(([m, ds]) =>
      `<div class="mono" style="font-size:.8rem">${K.escape(m)}: ${ds.map((d) => `${d.source} ${K.fmt(d.error_A, 1)} Å`).join(" · ")}</div>`).join("");
    return `<div class="flip" data-flipped="false" data-id="${K.escape(r.domain_id)}">
      <div class="flip-inner">
        <div class="flip-face">
          <div class="eyebrow mono">${K.escape(r.grp)}</div>
          <h3><a href="/k/${encodeURIComponent(r.domain_id)}">${K.escape(r.domain_id)}</a></h3>
          <p class="taxon">Degraded motif</p>
          <p class="mono" style="font-size:1.15rem;line-height:1.6">${K.escape(r.pseudo_evidence || "").split("; ").join("<br>")}</p>
          <button class="btn" type="button" data-flip aria-label="Show the geometry for ${K.escape(r.domain_id)}">Show geometry</button>
        </div>
        <div class="flip-face flip-back" inert>
          <h3>${K.escape(r.domain_id)}</h3>
          <ul style="list-style:none;padding:0;margin:0 0 8px;display:grid;gap:4px">${calls}</ul>
          ${errs || '<p class="muted">Geometry pending.</p>'}
          <div style="display:flex;gap:8px;flex-wrap:wrap;margin-top:10px">
            <button class="btn" type="button" data-flip aria-label="Show the sequence evidence for ${K.escape(r.domain_id)}">Show sequence</button>
            ${r.image_struct ? `<button class="btn" type="button" data-view="${K.escape(r.image_struct)}">View in 3D</button>` : ""}
          </div>
        </div>
      </div></div>`;
  }

  function render() {
    const f = document.getElementById("bs-filter").value;
    const show = rows.filter((r) => !f || (f === "restored" && (r.afdb_call === "restored" || r.esm_call === "restored")) ||
      (f === "faithful" && r.afdb_call !== "restored" && r.esm_call !== "restored" && (r.afdb_call || r.esm_call)) ||
      (f === "pdb" && r.pdb_call));
    document.getElementById("bs-count").textContent = `${show.length} of ${rows.length} pseudokinase domains`;
    const grid = document.getElementById("bs-grid");
    grid.innerHTML = show.map(card).join("");
    grid.querySelectorAll("[data-flip]").forEach((b) => b.addEventListener("click", () => {
      const f2 = b.closest(".flip");
      const on = f2.dataset.flipped !== "true";
      f2.dataset.flipped = String(on);
      f2.querySelector(".flip-face:not(.flip-back)").inert = on;
      f2.querySelector(".flip-back").inert = !on;
      (on ? f2.querySelector(".flip-back [data-flip]") : f2.querySelector(".flip-face [data-flip]")).focus();
    }));
    grid.querySelectorAll("[data-view]").forEach((b) => b.addEventListener("click", () => view(b.dataset.view, b.closest(".flip").dataset.id)));
  }

  async function view(path, id) {
    const dlg = document.getElementById("bs-dialog");
    document.getElementById("bs-dialog-title").textContent = `${id} · 3D view`;
    document.getElementById("bs-dialog-note").textContent = "Pocket-core frame; drag to rotate, scroll to zoom.";
    dlg.showModal();
    if (!viewer) viewer = await molstar.Viewer.create(document.getElementById("bs-molstar"), {
      layoutIsExpanded: false, layoutShowControls: false, layoutShowSequence: false, layoutShowLog: false,
      layoutShowLeftPanel: false, viewportShowExpand: false, viewportShowSelectionMode: false, viewportShowAnimation: false,
    });
    await viewer.plugin.clear();
    await viewer.loadStructureFromUrl(`/structures/${path}`, "mmcif", false);
  }

  document.getElementById("bs-close").addEventListener("click", () => document.getElementById("bs-dialog").close());
  document.getElementById("bs-filter").addEventListener("change", () => { K.state.set({ show: document.getElementById("bs-filter").value || null }); render(); });
  document.getElementById("bs-filter").value = K.state.get("show", "");
  K.json("/api/pseudo").then((r) => { rows = r; render(); });
})();
