/* Plate II: Specimen Card. */
(function () {
  "use strict";
  const K = window.KINFOLK;
  const C = K.colours;
  const GOBSMACKED_LINK_ENABLED = false; // see DECISIONS.md: GOBSMACKED does not yet read ?uniprot=&pocket=
  const domainId = document.getElementById("specimen").dataset.domain;
  const SRC_COLOUR = { pdb: 0x52625a, afdb: 0xc98f2b, esm2: 0xe07a55 };
  let data = null, viewer = null, ribbonCells = null, residueToPos = {};

  const devClass = (v) => (v === null || v === undefined ? "none" : v < 1 ? "ok" : v <= 2.5 ? "mid" : "high");
  const DEV_FILL = { ok: "#D9DFD0", mid: "#E6C58A", high: "#C4704F", none: "#FFFFFF" };

  // ------------------------------------------------------------ header bits
  function header() {
    document.getElementById("observer-note").textContent = data.observer_note;
    const fm = data.field_marks;
    const names = { b3_lys: "β3 Lys", ac_glu: "αC Glu", gatekeeper: "Gatekeeper", dfg_asp: "DFG Asp" };
    document.getElementById("field-marks").innerHTML = Object.entries(names).map(([k, label]) => {
      const f = fm[k] || {};
      return `<div><dt>${label}</dt><dd>${f.aa ? K.escape(f.aa) + f.resnum : "n/a"}</dd></div>`;
    }).join("");
    const a = data.agreement;
    const chip = (label, v, shape) =>
      `<span class="chip ${v === "agree" || v === "partial" || v === "disagree" ? v : ""}">${shape} ${label}: ${K.escape(K.agreementLabel[v] || v || "n/a")}</span>`;
    document.getElementById("agreement-chips").innerHTML =
      `<span class="chip">○ PDB states: ${a.pdb_states.length ? a.pdb_states.map(K.escape).join(", ") : "none"}</span>` +
      chip("AlphaFold DB", a.afdb, "■") + chip("ESMFold2", a.esm2, "▲");

    const gob = document.getElementById("gob-link");
    const residues = [fm.b3_lys, fm.ac_glu, fm.gatekeeper, fm.dfg_asp].filter((f) => f && f.resnum).map((f) => f.resnum);
    if (GOBSMACKED_LINK_ENABLED) {
      gob.href = `https://gobsmacked.mdeller.com/?uniprot=${encodeURIComponent(data.uniprot)}&pocket=${residues.join(",")}`;
      gob.removeAttribute("aria-disabled");
      document.getElementById("gob-wrap").removeAttribute("title");
    } else {
      gob.addEventListener("click", (e) => e.preventDefault());
    }

    const sel = document.getElementById("compare");
    const opts = [...(data.lookalikes.pocket || []).map((l) => ["pocket", l.domain]), ...(data.lookalikes.sequence || []).map((l) => ["sequence", l.domain])];
    const seen = new Set();
    opts.forEach(([kind, d]) => {
      if (seen.has(d)) return;
      seen.add(d);
      sel.insertAdjacentHTML("beforeend", `<option value="${K.escape(d)}">${K.escape(d)} (${kind} kin)</option>`);
    });
    sel.addEventListener("change", () => sel.value && (location.href = `/k/${encodeURIComponent(sel.value)}`));
  }

  // ------------------------------------------------------------ compass
  function compass() {
    const el = document.getElementById("compass");
    const W = 460, H = 360, m = { t: 14, r: 14, b: 40, l: 46 };
    const pts = [...data.compass.family.map(([d1, d2]) => ({ d1, d2 })), ...data.compass.pdb, ...data.compass.afdb, ...data.compass.esm]
      .filter((p) => p.d1 !== null && p.d2 !== null);
    const x = d3.scaleLinear().domain([Math.min(4, d3.min(pts, (p) => p.d1) || 4), Math.max(18, d3.max(pts, (p) => p.d1) || 18)]).nice().range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain([Math.min(4, d3.min(pts, (p) => p.d2) || 4), Math.max(18, d3.max(pts, (p) => p.d2) || 18)]).nice().range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("role", "img")
      .attr("aria-label", `Conformational compass for ${data.gene}: ${data.compass.pdb.length} PDB chains, AlphaFold DB model and ${data.compass.esm.length} ESMFold2 samples`);
    // Kincore regions
    const reg = svg.append("g").attr("fill", "#E9EDE2");
    reg.append("rect").attr("x", x(x.domain()[0])).attr("y", y(y.domain()[1])).attr("width", x(11) - x(x.domain()[0])).attr("height", y(11) - y(y.domain()[1]));
    reg.append("rect").attr("x", x(11)).attr("y", y(14)).attr("width", x(x.domain()[1]) - x(11)).attr("height", y(y.domain()[0]) - y(14)).attr("fill", "#F2E3DA");
    const lab = svg.append("g").attr("font-family", "Instrument Serif, serif").attr("font-style", "italic").attr("font-size", 15).attr("fill", "#52625A");
    lab.append("text").attr("x", x(x.domain()[0]) + 8).attr("y", y(y.domain()[1]) + 18).text("DFGin");
    lab.append("text").attr("x", x(x.domain()[0]) + 8).attr("y", y(y.domain()[0]) - 8).text("DFGinter");
    lab.append("text").attr("x", x(x.domain()[1]) - 8).attr("y", y(y.domain()[0]) - 8).attr("text-anchor", "end").text("DFGout");
    svg.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6)).call((g) => g.selectAll("text").attr("font-family", "JetBrains Mono, monospace"));
    svg.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6)).call((g) => g.selectAll("text").attr("font-family", "JetBrains Mono, monospace"));
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").attr("font-size", 12).attr("fill", C.ink).text("D1 (Å)");
    svg.append("text").attr("transform", `translate(12,${(H - m.b + m.t) / 2}) rotate(-90)`).attr("text-anchor", "middle").attr("font-size", 12).attr("fill", C.ink).text("D2 (Å)");
    svg.append("g").selectAll("circle").data(data.compass.family).enter().append("circle")
      .attr("cx", (d) => x(d[0])).attr("cy", (d) => y(d[1])).attr("r", 2).attr("fill", "#9AA39A").attr("opacity", 0.45);
    const tip = (sel, fmt) => sel.append("title").text(fmt);
    const gE = svg.append("g");
    gE.selectAll("path").data(data.compass.esm.filter((p) => p.d1 !== null)).enter().append("path")
      .attr("d", (p) => `M${x(p.d1)},${y(p.d2) - 6}l5.5,9.5h-11z`).attr("fill", C.esm).attr("stroke", C.ink).attr("stroke-width", 0.6)
      .attr("opacity", (p) => 0.25 + 0.75 * Math.max(0, Math.min(1, ((p.plddt || 50) - 50) / 45)))
      .call((s) => tip(s, (p) => `ESMFold2 sample ${p.i}: ${p.state}, pLDDT ${K.fmt(p.plddt, 1)}`));
    svg.append("g").selectAll("circle").data(data.compass.pdb.filter((p) => p.d1 !== null)).enter().append("circle")
      .attr("cx", (p) => x(p.d1)).attr("cy", (p) => y(p.d2)).attr("r", 5).attr("fill", "none").attr("stroke", C.ink).attr("stroke-width", 1.6)
      .call((s) => tip(s, (p) => `PDB ${p.pdb} ${p.chain}: ${p.state}${p.lig ? ", ligand " + p.lig : ""}`));
    svg.append("g").selectAll("rect").data(data.compass.afdb.filter((p) => p.d1 !== null)).enter().append("rect")
      .attr("x", (p) => x(p.d1) - 5).attr("y", (p) => y(p.d2) - 5).attr("width", 10).attr("height", 10)
      .attr("fill", C.afdb).attr("stroke", C.ink).attr("stroke-width", 1)
      .call((s) => tip(s, (p) => `AlphaFold DB: ${p.state}, pLDDT ${K.fmt(p.plddt, 1)}`));
  }

  // ------------------------------------------------------------ ribbon
  function ribbon() {
    const el = document.getElementById("ribbon");
    const rows = [["afdb", "AlphaFold DB", "■"], ["esm_best", "ESMFold2 best", "▲"], ["esm_mean", "ESMFold2 ensemble mean", "▲"]]
      .filter(([k]) => data.ribbon.rows[k]);
    if (!rows.length) {
      el.innerHTML = `<p class="muted">${data.is_dark ? "No PDB reference: deviations need an experimental structure." : "Predicted models pending."}</p>`;
      return;
    }
    const cw = 13, ch = 22, left = 170, top = 34;
    const W = left + 85 * cw + 10, H = top + rows.length * (ch + 4) + 40;
    const svg = d3.select(el).append("svg").attr("width", W).attr("height", H).attr("role", "img")
      .attr("aria-label", "KLIFS ribbon: deviation per pocket position; pale under 1 Å, amber 1 to 2.5 Å, rust over 2.5 Å");
    const pos = data.ribbon.positions;
    // region labels
    const regions = [];
    pos.forEach((p) => { const r = regions[regions.length - 1]; if (r && r.name === p.region) r.end = p.pos; else regions.push({ name: p.region, start: p.pos, end: p.pos }); });
    const rg = svg.append("g").attr("font-family", "JetBrains Mono, monospace").attr("font-size", 9).attr("fill", "#52625A");
    regions.forEach((r, i) => {
      const x0 = left + (r.start - 1) * cw, x1 = left + r.end * cw;
      rg.append("line").attr("x1", x0).attr("x2", x1 - 2).attr("y1", top - 6).attr("y2", top - 6).attr("stroke", i % 2 ? "#52625A" : "#1F3326");
      if (x1 - x0 >= 26) rg.append("text").attr("x", (x0 + x1) / 2).attr("y", top - 11).attr("text-anchor", "middle").text(r.name);
    });
    const cells = [];
    rows.forEach(([key, label, glyph], ri) => {
      const y0 = top + ri * (ch + 4);
      svg.append("text").attr("x", left - 8).attr("y", y0 + ch / 2 + 4).attr("text-anchor", "end").attr("font-size", 12).attr("fill", C.ink).text(`${glyph} ${label}`);
      const row = data.ribbon.rows[key];
      row.dev.forEach((v, i) => {
        const resnum = row.residue ? row.residue[i] : null;
        const r = svg.append("rect").attr("x", left + i * cw).attr("y", y0).attr("width", cw - 1).attr("height", ch)
          .attr("fill", DEV_FILL[devClass(v)]).attr("stroke", "#C9D0C1").attr("stroke-width", 0.5)
          .attr("data-pos", i + 1).attr("data-row", key);
        r.append("title").text(`KLIFS ${i + 1} (${pos[i].region}${pos[i].aa ? ", " + pos[i].aa : ""}): ${v === null ? "no data" : v.toFixed(2) + " Å"}`);
        r.on("mouseenter", () => highlight(i + 1, key));
        cells.push(r);
        if (key !== "esm_mean" && resnum !== null && resnum !== undefined) (residueToPos[key] = residueToPos[key] || {})[resnum] = i + 1;
      });
    });
    const ly = top + rows.length * (ch + 4) + 16;
    [["ok", "< 1 Å"], ["mid", "1 to 2.5 Å"], ["high", "> 2.5 Å"], ["none", "missing"]].forEach(([k, t], i) => {
      svg.append("rect").attr("x", left + i * 110).attr("y", ly - 10).attr("width", 14).attr("height", 14).attr("fill", DEV_FILL[k]).attr("stroke", "#C9D0C1");
      svg.append("text").attr("x", left + i * 110 + 20).attr("y", ly + 1).attr("font-size", 12).attr("fill", C.ink).text(t);
    });
    ribbonCells = svg;
  }

  function markCell(pos) {
    if (!ribbonCells) return;
    ribbonCells.selectAll("rect[data-pos]").attr("stroke", "#C9D0C1").attr("stroke-width", 0.5);
    ribbonCells.selectAll(`rect[data-pos="${pos}"]`).attr("stroke", C.disagree).attr("stroke-width", 2.6).raise();
  }

  // ------------------------------------------------------------ Mol*
  function findNs(obj, key, depth = 0) {
    if (!obj || depth > 3 || typeof obj !== "object") return null;
    if (obj[key]) return obj[key];
    for (const k of Object.keys(obj)) {
      const r = findNs(obj[k], key, depth + 1);
      if (r) return r;
    }
    return null;
  }

  function pdbLabel(id) {
    const p = (data.compass.pdb || []).find((x) => x.id === id);
    return p ? `${p.pdb} ${p.chain}` : "n/a";
  }

  function structuresToLoad() {
    const s = data.structures || {};
    const ids = Object.keys(s);
    const pdb = ids.filter((i) => i.startsWith("pdb:"));
    const af = ids.filter((i) => i.startsWith("afdb:"));
    const esm = ids.filter((i) => i.startsWith("esm2:"));
    return { pdb: pdb[0], afdb: af[0], esm };
  }

  async function initMolstar() {
    const host = document.getElementById("molstar");
    const set = structuresToLoad();
    const esmSel = document.getElementById("mol-esm");
    set.esm.forEach((id) => esmSel.insertAdjacentHTML("beforeend", `<option value="${K.escape(id)}">${K.escape(id.split(":").slice(1).join(" · "))}</option>`));
    if (!set.esm.length) esmSel.insertAdjacentHTML("beforeend", `<option value="">No ESMFold2 samples yet</option>`);
    if (!window.molstar || (!set.pdb && !set.afdb && !set.esm.length)) {
      host.innerHTML = `<p class="muted" style="padding:16px">No structures bundled for this domain.</p>`;
      return;
    }
    viewer = await molstar.Viewer.create(host, {
      layoutIsExpanded: false, layoutShowControls: false, layoutShowSequence: false, layoutShowLog: false,
      layoutShowLeftPanel: false, viewportShowExpand: false, viewportShowSelectionMode: false,
      viewportShowAnimation: false, viewportShowTrajectoryControls: false, collapseLeftPanel: true,
      pdbProvider: "rcsb", emdbProvider: "rcsb",
    });
    await loadSet(set.pdb, set.afdb, set.esm[0]);
    esmSel.addEventListener("change", () => loadSet(set.pdb, set.afdb, esmSel.value));
    subscribeClicks();
    document.getElementById("mol-legend").innerHTML =
      `<span>○ grey: PDB ${K.escape(pdbLabel(set.pdb))}</span><span>■ amber: AlphaFold DB</span><span>▲ coral: ESMFold2</span>`;
  }

  async function loadSet(pdb, afdb, esm) {
    await viewer.plugin.clear();
    const load = (id, src) => id && viewer.loadStructureFromUrl(data.structures[id], "mmcif", false, {
      representationParams: { theme: { globalName: "uniform", globalColorParams: { value: SRC_COLOUR[src] } } },
      label: id,
    });
    await load(pdb, "pdb");
    await load(afdb, "afdb");
    await load(esm, "esm2");
  }

  function highlight(pos, rowKey) {
    markCell(pos);
    if (!viewer) return;
    const key = rowKey === "esm_mean" ? "esm_best" : rowKey;
    const row = data.ribbon.rows[key] || data.ribbon.rows.afdb;
    const res = row && row.residue ? row.residue[pos - 1] : null;
    if (res === null || res === undefined) return;
    try {
      viewer.structureInteractivity({ elements: { auth_seq_id: res }, action: "highlight" });
    } catch (e) { /* viewer not ready */ }
  }

  function subscribeClicks() {
    const SE = findNs(molstar.lib, "StructureElement");
    const SP = findNs(molstar.lib, "StructureProperties");
    if (!SE || !SP) return;
    viewer.plugin.behaviors.interaction.click.subscribe((ev) => {
      const loci = ev.current && ev.current.loci;
      if (!loci || !SE.Loci.is(loci) || SE.Loci.isEmpty(loci)) return;
      const loc = SE.Location.create();
      SE.Loci.getFirstLocation(loci, loc);
      const resnum = SP.residue.auth_seq_id(loc);
      const label = (loc.structure && loc.structure.label) || "";
      const key = label.startsWith("afdb") ? "afdb" : "esm_best";
      const pos = (residueToPos[key] || residueToPos.afdb || {})[resnum];
      if (pos) {
        markCell(pos);
        document.getElementById("ribbon").scrollTo({ left: Math.max(0, 170 + (pos - 8) * 13), behavior: "smooth" });
      }
    });
  }

  // ------------------------------------------------------------ cards
  function cards() {
    const rooms = document.getElementById("card-rooms");
    rooms.innerHTML = data.sites.length
      ? `<ul style="padding-left:18px;margin:0">${data.sites.slice(0, 6).map((s) =>
          `<li><span class="mono">${K.escape(s.id.split(":")[1])}</span> ${K.escape(s.subpocket)}: recovered in ${Math.round((s.recovery || 0) * 100)}% of samples${s.candidate ? " · <b>candidate cryptic</b>" : s.known ? " · known" : ""}</li>`).join("")}</ul>
         <p><a href="/rooms?k=${encodeURIComponent(data.domain_id)}">Open in Hidden Rooms</a></p>`
      : `<p class="muted">No pocket sites detected yet.</p>`;
    const look = document.getElementById("card-look");
    const pl = data.lookalikes.pocket || [], sl = data.lookalikes.sequence || [];
    look.innerHTML = `<div class="grid-2" style="gap:10px">
      <div><div class="ctl-label">Pocket kin</div><ol style="padding-left:18px;margin:0">${pl.map((l) => `<li><a href="/k/${encodeURIComponent(l.domain)}">${K.escape(l.domain)}</a> <span class="mono muted">${K.fmt(l.sim, 2)}</span></li>`).join("") || "<li>n/a</li>"}</ol></div>
      <div><div class="ctl-label">Sequence kin</div><ol style="padding-left:18px;margin:0">${sl.map((l) => `<li><a href="/k/${encodeURIComponent(l.domain)}">${K.escape(l.domain)}</a> <span class="mono muted">${K.fmt(l.dist, 2)}</span></li>`).join("") || "<li>n/a</li>"}</ol></div></div>`;
    const r = data.range;
    document.getElementById("card-range").innerHTML = `<dl class="field-marks" style="grid-template-columns:1fr 1fr">
      <div><dt>PDB chains</dt><dd>${r.n_chains}</dd></div><div><dt>Blind (ESMFold2)</dt><dd>${r.n_blind_esm}</dd></div>
      <div><dt>Blind + novel</dt><dd>${r.n_novel_esm}</dd></div><div><dt>First release</dt><dd>${K.escape(r.first_release || "n/a")}</dd></div></dl>
      <p style="margin:10px 0 0;font-size:.88rem">${(r.ligands || []).slice(0, 8).map((l) => `<span class="mono">${K.escape(l.ccd)}</span> <span class="muted">(type ${K.escape(l.type || "?")}, ${l.n})</span>`).join(" · ") || "No ligands recorded."}</p>`;
  }

  K.json(`/static/data/kinase/${encodeURIComponent(domainId)}.json`).then((d) => {
    data = d;
    header();
    compass();
    ribbon();
    cards();
    initMolstar().catch((e) => {
      document.getElementById("molstar").innerHTML = `<p class="muted" style="padding:16px">The 3D viewer could not start (${K.escape(e.message)}).</p>`;
    });
  });
})();
