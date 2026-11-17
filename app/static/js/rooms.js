/* Plate IV: Hidden Rooms (site atlas, site viewer, sunburst). */
(function () {
  "use strict";
  const K = window.KINFOLK;
  const SUB_COLOURS = { front: "#4F6C82", gate: "#C98F2B", "back pocket BP-I": "#9A5C7A", "back pocket BP-II": "#A9482A",
    "outside KLIFS (allosteric)": "#3F6B52" };
  let sites = [], sortKey = "rank_score", sortDir = -1, viewer = null, current = null, domainSites = [];

  function filtered() {
    const g = document.getElementById("f-group").value, s = document.getElementById("f-sub").value;
    const cand = document.getElementById("f-cand").value === "1";
    return sites.filter((r) => (!g || r.grp === g) && (!s || r.subpocket === s) && (!cand || r.candidate))
      .sort((a, b) => sortDir * (((a[sortKey] ?? -1) > (b[sortKey] ?? -1)) - ((a[sortKey] ?? -1) < (b[sortKey] ?? -1))));
  }

  function table() {
    const rows = filtered();
    document.querySelector("#sites tbody").innerHTML = rows.slice(0, 400).map((r) => `<tr>
      <td><a href="/k/${encodeURIComponent(r.domain_id)}">${K.escape(r.domain_id)}</a></td><td>${K.escape(r.grp)}</td>
      <td>${K.escape(r.subpocket)}</td><td class="num">${K.fmt(r.recovery, 2)}</td><td class="num">${K.fmt(r.p2rank_max, 2)}</td>
      <td class="num">${r.vol_max ? Math.round(r.vol_max) : "n/a"}</td><td class="num">${K.fmt(r.rank_score, 3)}</td>
      <td><button class="btn" type="button" data-site="${K.escape(r.site_id)}">View</button></td></tr>`).join("") ||
      `<tr><td colspan="8" class="muted">No sites match.</td></tr>`;
    document.querySelectorAll("#sites tbody button[data-site]").forEach((b) => b.addEventListener("click", () => show(b.dataset.site)));
  }

  function sunburst() {
    const cand = sites.filter((r) => r.candidate);
    const el = document.getElementById("sunburst");
    if (!cand.length) { el.innerHTML = `<p class="muted">No candidate sites yet.</p>`; return; }
    const nest = { name: "sites" };
    const groups = d3.rollup(cand, (v) => v.length, (r) => r.subpocket, (r) => r.grp);
    nest.children = [...groups].map(([sp, m]) => ({ name: sp, children: [...m].map(([g, n]) => ({ name: g, value: n })) }));
    const R = 170;
    const root = d3.partition().size([2 * Math.PI, R])(d3.hierarchy(nest).sum((d) => d.value || 0));
    const arc = d3.arc().startAngle((d) => d.x0).endAngle((d) => d.x1).innerRadius((d) => d.y0).outerRadius((d) => d.y1 - 1).padAngle(0.004);
    const svg = d3.select(el).append("svg").attr("viewBox", `${-R - 10} ${-R - 10} ${2 * R + 20} ${2 * R + 20}`);
    svg.append("g").selectAll("path").data(root.descendants().filter((d) => d.depth)).enter().append("path").attr("d", arc)
      .attr("fill", (d) => SUB_COLOURS[d.depth === 1 ? d.data.name : d.parent.data.name] || "#9AA39A")
      .attr("opacity", (d) => (d.depth === 1 ? 1 : 0.55)).attr("stroke", "#F8F7F1")
      .append("title").text((d) => `${d.depth === 1 ? d.data.name : d.parent.data.name + " › " + d.data.name}: ${d.value} sites`);
    svg.append("g").attr("font-size", 11).attr("text-anchor", "middle").attr("pointer-events", "none")
      .selectAll("text").data(root.descendants().filter((d) => d.depth === 1 && d.x1 - d.x0 > 0.25)).enter().append("text")
      .attr("transform", (d) => { const a = ((d.x0 + d.x1) / 2) * 180 / Math.PI; const r = (d.y0 + d.y1) / 2; return `rotate(${a - 90}) translate(${r},0) rotate(${a < 180 ? 0 : 180})`; })
      .attr("fill", "#fff").text((d) => d.data.name.replace("back pocket ", "").replace(" (allosteric)", ""));
    svg.append("text").attr("text-anchor", "middle").attr("dy", "0.35em").attr("font-family", "Instrument Serif, serif").attr("font-size", 22).text(`${cand.length}`);
  }

  // ------------------------------------------------------------ viewer
  function findNs(obj, key, depth = 0) {
    if (!obj || depth > 3 || typeof obj !== "object") return null;
    if (obj[key]) return obj[key];
    for (const k of Object.keys(obj)) { const r = findNs(obj[k], key, depth + 1); if (r) return r; }
    return null;
  }

  async function show(siteId) {
    const site = sites.find((s) => s.site_id === siteId);
    if (!site) return;
    K.state.set({ site: siteId });
    const spec = await K.json(`/static/data/kinase/${encodeURIComponent(site.domain_id)}.json`);
    const esm = Object.keys(spec.structures || {}).find((k) => k.startsWith("esm2:"));
    const url = esm ? spec.structures[esm] : Object.values(spec.structures || {})[0];
    document.getElementById("viewer-label").textContent = `${site.domain_id} · ${esm || "no ESMFold2 model"}`;
    domainSites = (await K.json(`/api/sites/${encodeURIComponent(site.domain_id)}`)).map((s) => ({ ...s, res: JSON.parse(s.residues_json) }));
    current = { site, url, esm, construct: spec.dom_start };
    if (!viewer) {
      viewer = await molstar.Viewer.create(document.getElementById("site-molstar"), {
        layoutIsExpanded: false, layoutShowControls: false, layoutShowSequence: false, layoutShowLog: false,
        layoutShowLeftPanel: false, viewportShowExpand: false, viewportShowSelectionMode: false, viewportShowAnimation: false,
      });
    }
    await draw();
  }

  async function draw() {
    if (!viewer || !current || !current.url) return;
    const thr = +document.getElementById("thresh").value;
    await viewer.plugin.clear();
    await viewer.loadStructureFromUrl(current.url, "mmcif", false, {
      representationParams: { theme: { globalName: "uniform", globalColorParams: { value: 0xd9dfd0 } } },
    });
    const MS = findNs(molstar.lib, "MolScriptBuilder");
    const plugin = viewer.plugin;
    const struct = plugin.managers.structure.hierarchy.current.structures[0];
    const show = domainSites.filter((s) => (s.recovery || 0) >= thr).slice(0, 6);
    const legend = [];
    const palette = [0xe07a55, 0x3f6b52, 0xc98f2b, 0x4f6c82, 0x9a5c7a, 0xa9482a];
    // ESMFold2 models are numbered from the construct start; site residues are UniProt numbers.
    const offset = current.esm ? (await constructOffset(current.site.domain_id, current.esm)) : 0;
    for (const [i, s] of show.entries()) {
      const nums = s.res.residues.map((u) => u - offset);
      if (!MS || !struct) break;
      const expr = MS.struct.generator.atomGroups({
        "residue-test": MS.core.set.has([MS.set(...nums), MS.struct.atomProperty.macromolecular.auth_seq_id()]),
      });
      try {
        const comp = await plugin.builders.structure.tryCreateComponentFromExpression(struct.cell, expr, `site-${i}`, { label: s.site_id });
        if (comp) await plugin.builders.structure.representation.addRepresentation(comp, {
          type: "molecular-surface", typeParams: { alpha: s.site_id === current.site.site_id ? 0.7 : 0.35 },
          color: "uniform", colorParams: { value: palette[i % palette.length] },
        });
      } catch (e) { /* representation unavailable: keep cartoon */ }
      legend.push(`<span><svg width="12" height="12"><rect width="12" height="12" fill="#${palette[i % palette.length].toString(16)}"/></svg>${K.escape(s.site_id.split(":")[1])} · ${K.escape(s.subpocket)} · ${Math.round((s.recovery || 0) * 100)}%</span>`);
    }
    document.getElementById("site-legend").innerHTML = legend.join("") || "<span>No sites above this threshold.</span>";
  }

  const offsets = new Map();
  async function constructOffset(domainId, esmId) {
    if (offsets.has(esmId)) return offsets.get(esmId);
    const rows = await K.json(`/api/summary/construct_starts`).catch(() => ({}));
    const v = (rows[domainId] || 1) - 1;
    offsets.set(esmId, v);
    return v;
  }

  // ------------------------------------------------------------ wiring
  document.querySelectorAll("#sites thead button").forEach((b) => b.addEventListener("click", () => {
    const k = b.dataset.k;
    sortDir = sortKey === k ? -sortDir : -1;
    sortKey = k;
    document.querySelectorAll("#sites thead button").forEach((x) => x.closest("th").removeAttribute("aria-sort"));
    b.closest("th").setAttribute("aria-sort", sortDir < 0 ? "descending" : "ascending");
    table();
  }));
  ["f-group", "f-sub", "f-cand"].forEach((id) => document.getElementById(id).addEventListener("change", table));
  const th = document.getElementById("thresh");
  th.addEventListener("input", () => { document.getElementById("thresh-val").textContent = (+th.value).toFixed(2); });
  th.addEventListener("change", draw);

  K.json("/api/sites").then((rows) => {
    sites = rows;
    const groups = [...new Set(rows.map((r) => r.grp))].sort(), subs = [...new Set(rows.map((r) => r.subpocket))].sort();
    document.getElementById("f-group").insertAdjacentHTML("beforeend", groups.map((g) => `<option>${K.escape(g)}</option>`).join(""));
    document.getElementById("f-sub").insertAdjacentHTML("beforeend", subs.map((g) => `<option>${K.escape(g)}</option>`).join(""));
    table();
    sunburst();
    const pre = K.state.get("site", null) || (K.state.get("k", null) && (rows.find((r) => r.domain_id === K.state.get("k")) || {}).site_id);
    if (pre) show(pre);
  });
})();
