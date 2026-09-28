/* Plate VI: dark kinome table and agreement scatter. */
(function () {
  "use strict";
  const K = window.KINFOLK;
  let rows = [], sortKey = "esm_afdb_rmsd", dir = 1;

  function filtered() {
    const g = document.getElementById("dk-group").value, q = document.getElementById("dk-q").value.trim().toLowerCase();
    return rows.filter((r) => (!g || r.grp === g) && (!q || `${r.domain_id} ${r.family}`.toLowerCase().includes(q)))
      .sort((a, b) => {
        const x = a[sortKey], y = b[sortKey];
        if (x === null || x === undefined) return 1;
        if (y === null || y === undefined) return -1;
        return dir * (x > y ? 1 : x < y ? -1 : 0);
      });
  }

  function table() {
    const f = filtered();
    document.getElementById("dk-count").textContent = `${f.length} of ${rows.length} dark domains`;
    document.querySelector("#dk-table tbody").innerHTML = f.map((r) => `<tr>
      <td><a href="/k/${encodeURIComponent(r.domain_id)}">${K.escape(r.domain_id)}</a></td><td>${K.escape(r.grp)}</td><td>${K.escape(r.family)}</td>
      <td class="num">${K.fmt(r.esm_afdb_rmsd, 2)}</td><td class="num">${K.fmt(r.esm_entropy, 2)}</td>
      <td class="num">${K.fmt(r.afdb_plddt, 1)}</td><td class="num">${K.fmt(r.esm_plddt, 1)}</td>
      <td>${r.top_site ? `<a href="/rooms?site=${encodeURIComponent(r.top_site)}">${K.escape(r.top_site.split(":")[1])}</a>` : "n/a"}</td>
      <td class="num">${r.n_candidates}</td></tr>`).join("");
    scatter(f);
  }

  function scatter(f) {
    const el = document.getElementById("dk-scatter");
    el.innerHTML = "";
    const pts = f.filter((r) => r.esm_afdb_rmsd !== null && r.esm_entropy !== null);
    if (!pts.length) { el.innerHTML = `<p class="muted">ESMFold2 ensembles for the dark kinome are still folding.</p>`; return; }
    const W = 460, H = 320, m = { t: 12, r: 16, b: 42, l: 50 };
    const x = d3.scaleLinear().domain([0, d3.max(pts, (p) => p.esm_afdb_rmsd) || 1]).nice().range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain([0, Math.max(1, d3.max(pts, (p) => p.esm_entropy))]).nice().range([H - m.b, m.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
    svg.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6));
    svg.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
    svg.append("text").attr("x", (m.l + W) / 2).attr("y", H - 6).attr("text-anchor", "middle").attr("font-size", 12).text("ESMFold2 vs AlphaFold DB pocket RMSD (Å)");
    svg.append("text").attr("transform", `translate(12,${H / 2}) rotate(-90)`).attr("text-anchor", "middle").attr("font-size", 12).text("Ensemble state entropy (bits)");
    svg.append("g").selectAll("a").data(pts).enter().append("a").attr("href", (p) => `/k/${encodeURIComponent(p.domain_id)}`)
      .append("path").attr("d", (p) => `M${x(p.esm_afdb_rmsd)},${y(p.esm_entropy) - 5 - p.n_candidates}l${5 + p.n_candidates},${9 + 2 * p.n_candidates}h${-10 - 2 * p.n_candidates}z`)
      .attr("fill", "#E07A55").attr("stroke", "#1F3326").attr("stroke-width", 0.6).attr("opacity", 0.8)
      .append("title").text((p) => `${p.domain_id}: ${K.fmt(p.esm_afdb_rmsd, 2)} Å, ${K.fmt(p.esm_entropy, 2)} bits, ${p.n_candidates} candidate sites`);
  }

  document.querySelectorAll("#dk-table thead button").forEach((b) => b.addEventListener("click", () => {
    dir = sortKey === b.dataset.k ? -dir : 1;
    sortKey = b.dataset.k;
    document.querySelectorAll("#dk-table thead button").forEach((x) => x.removeAttribute("aria-sort"));
    b.setAttribute("aria-sort", dir > 0 ? "ascending" : "descending");
    table();
  }));
  ["dk-group", "dk-q"].forEach((id) => document.getElementById(id).addEventListener("input", table));
  K.json("/api/dark").then((r) => {
    rows = r;
    const groups = [...new Set(r.map((x) => x.grp))].sort();
    document.getElementById("dk-group").insertAdjacentHTML("beforeend", groups.map((g) => `<option>${K.escape(g)}</option>`).join(""));
    table();
  });
})();
