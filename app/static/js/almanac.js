/* Plate VII: The Almanac (Q1 to Q6 figures from live summaries). */
(function () {
  "use strict";
  const K = window.KINFOLK;
  const C = K.colours;
  const D = JSON.parse(document.getElementById("alm-data").textContent || "{}");
  const F = D.findings || {};
  const STATE_COLOURS = ["#3F6B52", "#6E8F5B", "#8A7A3A", "#4F6C82", "#7C6A92", "#9A5C7A", "#C98F2B", "#A9482A", "#9AA39A"];
  const pending = (el, what) => (el.innerHTML = `<p class="muted">${what} pending: the fold queue is still running.</p>`);

  document.querySelectorAll("[data-finding]").forEach((p) => {
    const f = F[p.dataset.finding];
    p.textContent = f && f.summary ? f.summary : "Result pending.";
  });

  // ---------------------------------------------------------------- Q1
  K.json("/api/summary/motif_fidelity").then((m) => {
    const el = document.getElementById("fig-q1");
    if (!m || !m.afdb || !Object.keys(m.afdb).length) return pending(el, "Motif deviations");
    const motifs = ["g.l", "aC", "b.l", "GK", "hinge", "c.l", "xDFG", "a.l"];
    const labels = { "g.l": "G-loop", aC: "αC", "b.l": "β-loop", GK: "gatekeeper", hinge: "hinge", "c.l": "catalytic loop", xDFG: "xDFG", "a.l": "activation loop" };
    const W = 640, H = 300, mg = { t: 16, r: 16, b: 56, l: 50 };
    const max = d3.max(["afdb", "esm2_best"].flatMap((s) => motifs.map((mo) => (m[s] && m[s][mo] ? m[s][mo].q3 : 0)))) || 1;
    const x0 = d3.scaleBand().domain(motifs).range([mg.l, W - mg.r]).padding(0.25);
    const x1 = d3.scaleBand().domain(["afdb", "esm2_best"]).range([0, x0.bandwidth()]).padding(0.15);
    const y = d3.scaleLinear().domain([0, max]).nice().range([H - mg.b, mg.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
    svg.append("g").attr("transform", `translate(0,${H - mg.b})`).call(d3.axisBottom(x0).tickFormat((d) => labels[d]))
      .selectAll("text").attr("transform", "rotate(-25)").attr("text-anchor", "end");
    svg.append("g").attr("transform", `translate(${mg.l},0)`).call(d3.axisLeft(y).ticks(5));
    svg.append("text").attr("transform", `translate(12,${H / 2}) rotate(-90)`).attr("text-anchor", "middle").attr("font-size", 12).text("Cα RMSD to nearest PDB state (Å)");
    motifs.forEach((mo) => ["afdb", "esm2_best"].forEach((s) => {
      const v = m[s] && m[s][mo];
      if (!v) return;
      const cx = x0(mo) + x1(s) + x1.bandwidth() / 2;
      svg.append("line").attr("x1", cx).attr("x2", cx).attr("y1", y(v.q1)).attr("y2", y(v.q3)).attr("stroke", C.ink);
      if (s === "afdb") svg.append("rect").attr("x", cx - 5).attr("y", y(v.median) - 5).attr("width", 10).attr("height", 10).attr("fill", C.afdb).attr("stroke", C.ink)
        .append("title").text(`AlphaFold DB ${labels[mo]}: median ${v.median} Å (n ${v.n})`);
      else svg.append("path").attr("d", `M${cx},${y(v.median) - 7}l6,11h-12z`).attr("fill", C.esm).attr("stroke", C.ink)
        .append("title").text(`ESMFold2 best ${labels[mo]}: median ${v.median} Å (n ${v.n})`);
    }));
  }).catch(() => pending(document.getElementById("fig-q1"), "Motif deviations"));

  // ---------------------------------------------------------------- Q2 (PCI morph)
  (function q2() {
    const el = document.getElementById("fig-q2");
    const pci = D.pci || {};
    if (!pci.distributions) return pending(el, "The Prior Collapse Index");
    const srcs = [["pdb", "PDB"], ["afdb", "AlphaFold DB"], ["esm_single", "ESMFold2 best sample"], ["esm_ensemble", "ESMFold2 ensemble"]];
    const states = pci.states;
    const W = 640, H = 300, mg = { t: 16, r: 16, b: 60, l: 50 };
    const x = d3.scaleBand().domain(states).range([mg.l, W - mg.r]).padding(0.2);
    const y = d3.scaleLinear().domain([0, d3.max(srcs.flatMap(([s]) => pci.distributions[s]))]).nice().range([H - mg.b, mg.t]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
    svg.append("g").attr("transform", `translate(0,${H - mg.b})`).call(d3.axisBottom(x)).selectAll("text").attr("transform", "rotate(-30)").attr("text-anchor", "end");
    svg.append("g").attr("transform", `translate(${mg.l},0)`).call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".0%")));
    const ghost = svg.append("g");
    ghost.selectAll("rect").data(pci.distributions.pdb).enter().append("rect").attr("x", (d, i) => x(states[i])).attr("width", x.bandwidth())
      .attr("y", (d) => y(d)).attr("height", (d) => y(0) - y(d)).attr("fill", "none").attr("stroke", C.ink).attr("stroke-dasharray", "3 2");
    const bars = svg.append("g").selectAll("rect").data(pci.distributions.pdb).enter().append("rect")
      .attr("x", (d, i) => x(states[i])).attr("width", x.bandwidth()).attr("fill", (d, i) => STATE_COLOURS[i]).attr("opacity", 0.85)
      .attr("y", (d) => y(d)).attr("height", (d) => y(0) - y(d));
    const slider = document.getElementById("pci-slider");
    const update = () => {
      const [s, label] = srcs[+slider.value];
      document.getElementById("pci-label").textContent = label;
      slider.setAttribute("aria-valuetext", label);
      bars.data(pci.distributions[s]).transition().duration(600).attr("y", (d) => y(d)).attr("height", (d) => y(0) - y(d));
    };
    slider.addEventListener("input", update);
    const k = pci.kinome;
    document.getElementById("pci-table").innerHTML = `<div class="table-wrap"><table><caption class="sr-only">Prior Collapse Index</caption>
      <thead><tr><th scope="col">Comparison with the PDB</th><th scope="col" class="num">PCI (JSD, bits)</th><th scope="col" class="num">95% CI</th></tr></thead><tbody>
      ${[["afdb", "AlphaFold DB"], ["esm_single", "ESMFold2 best sample"], ["esm_ensemble", "ESMFold2 ensemble"]].map(([s, l]) =>
        `<tr><td>${l}</td><td class="num">${k[s].jsd.toFixed(3)}</td><td class="num">${k[s].ci95[0].toFixed(3)} to ${k[s].ci95[1].toFixed(3)}</td></tr>`).join("")}
      </tbody></table></div><p class="muted" style="font-size:.85rem">${pci.n_kinases} kinases with PDB, AlphaFold DB and ESMFold2 data; dashed outline shows the PDB distribution.</p>`;
  })();

  // ---------------------------------------------------------------- Q3
  (function q3() {
    const el = document.getElementById("fig-q3");
    const f = F.q3 && F.q3.blind;
    if (!f) return pending(el, "The blind-set comparison");
    const cats = Object.keys(f);
    const W = 640, H = 260, mg = { t: 16, r: 16, b: 50, l: 180 };
    const y = d3.scaleBand().domain(cats).range([mg.t, H - mg.b]).padding(0.3);
    const x = d3.scaleLinear().domain([0, d3.max(cats.flatMap((c) => ["afdb", "esm2"].map((s) => (f[c][s] ? f[c][s].q3 : 0))))]).nice().range([mg.l, W - mg.r]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
    svg.append("g").attr("transform", `translate(0,${H - mg.b})`).call(d3.axisBottom(x).ticks(6));
    svg.append("g").attr("transform", `translate(${mg.l},0)`).call(d3.axisLeft(y));
    svg.append("text").attr("x", (mg.l + W) / 2).attr("y", H - 8).attr("text-anchor", "middle").attr("font-size", 12).text("Pocket Cα RMSD to the blind structure (Å)");
    cats.forEach((c) => ["afdb", "esm2"].forEach((s, i) => {
      const v = f[c][s];
      if (!v) return;
      const cy = y(c) + y.bandwidth() * (i ? 0.72 : 0.28);
      svg.append("line").attr("x1", x(v.q1)).attr("x2", x(v.q3)).attr("y1", cy).attr("y2", cy).attr("stroke", C.ink);
      if (s === "afdb") svg.append("rect").attr("x", x(v.median) - 5).attr("y", cy - 5).attr("width", 10).attr("height", 10).attr("fill", C.afdb).attr("stroke", C.ink);
      else svg.append("path").attr("d", `M${x(v.median)},${cy - 7}l6,11h-12z`).attr("fill", C.esm).attr("stroke", C.ink);
    }));
  })();

  // ---------------------------------------------------------------- Q4
  K.json("/api/summary/controls").then((ctrl) => {
    const el = document.getElementById("fig-q4");
    if (!ctrl || !ctrl.length) return pending(el, "Positive controls");
    el.innerHTML = `<div class="table-wrap"><table><caption class="sr-only">Positive controls</caption><thead><tr><th scope="col">Control</th><th scope="col">Site</th>
      <th scope="col">▲ ESMFold2</th><th scope="col" class="num">Recovery</th><th scope="col">■ AlphaFold DB</th></tr></thead><tbody>
      ${ctrl.map((c) => `<tr><td>${K.escape(c.gene)}</td><td>${K.escape(c.site_name)}</td>
        <td>${c.recovered === null || c.recovered === undefined ? "pending" : c.recovered ? '<span class="chip agree">recovered</span>' : '<span class="chip disagree">not recovered</span>'}</td>
        <td class="num">${c.score === null || c.score === undefined ? "n/a" : Math.round(c.score * 100) + "%"}</td>
        <td>${c.afdb_recovered ? '<span class="chip agree">recovered</span>' : c.afdb_jaccard !== undefined ? `Jaccard ${K.fmt(c.afdb_jaccard, 2)}` : "n/a"}</td></tr>`).join("")}
      </tbody></table></div>`;
  }).catch(() => pending(document.getElementById("fig-q4"), "Positive controls"));

  // ---------------------------------------------------------------- Q5
  (function () {
    const q = D.q5 || {};
    const el = document.getElementById("fig-q5");
    if (!el) return;
    if (!q.rho) return pending(el, "Selectivity correlations");
    const srcs = [["pdb", "PDB", "circle"], ["afdb", "AlphaFold DB", "square"], ["esm2", "ESMFold2", "triangle"]]
      .filter(([s]) => q.rho[s] !== undefined);
    const W = 560, H = 60 + srcs.length * 54, mg = { l: 140, r: 40, t: 20, b: 40 };
    const lo = Math.min(0, ...srcs.map(([s]) => (q.rho_ci95[s] || [0])[0]));
    const hi = Math.max(0.1, ...srcs.map(([s]) => (q.rho_ci95[s] || [0, 0])[1]));
    const x = d3.scaleLinear().domain([lo, hi]).nice().range([mg.l, W - mg.r]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
    svg.append("g").attr("transform", `translate(0,${H - mg.b})`).call(d3.axisBottom(x).ticks(5));
    svg.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", mg.t).attr("y2", H - mg.b).attr("stroke", C.rule);
    svg.append("text").attr("x", (mg.l + W - mg.r) / 2).attr("y", H - 6).attr("text-anchor", "middle")
      .attr("font-size", 12).text("Spearman \u03c1 with Davis co-inhibition");
    srcs.forEach(([s, label, shape], i) => {
      const y = mg.t + 24 + i * 54, [a, b] = q.rho_ci95[s];
      svg.append("text").attr("x", mg.l - 10).attr("y", y + 4).attr("text-anchor", "end").attr("font-size", 13).text(label);
      svg.append("line").attr("x1", x(a)).attr("x2", x(b)).attr("y1", y).attr("y2", y).attr("stroke", C.ink).attr("stroke-width", 2);
      const cx = x(q.rho[s]);
      const fill = s === "pdb" ? "#fff" : s === "afdb" ? C.afdb : C.esm;
      if (shape === "circle") svg.append("circle").attr("cx", cx).attr("cy", y).attr("r", 7).attr("fill", fill).attr("stroke", C.ink).attr("stroke-width", 2);
      else if (shape === "square") svg.append("rect").attr("x", cx - 7).attr("y", y - 7).attr("width", 14).attr("height", 14).attr("fill", fill).attr("stroke", C.ink);
      else svg.append("path").attr("d", `M${cx},${y - 8}l7,12h-14z`).attr("fill", fill).attr("stroke", C.ink);
      svg.append("text").attr("x", x(b) + 8).attr("y", y + 4).attr("font-family", "JetBrains Mono, monospace")
        .attr("font-size", 11).text(q.rho[s].toFixed(2));
    });
    const f = q.fraction_of_pdb || {};
    const note = document.getElementById("q5-note");
    if (note) {
      note.textContent = `${q.n_kinases} kinases, ${q.n_pairs} pairs. Predicted-structure correlation as a fraction `
        + "of the PDB correlation: "
        + Object.entries(f).map(([s, v]) =>
            `${s === "esm2" ? "ESMFold2" : "AlphaFold DB"} ${v.toFixed(2)} `
            + `(95% CI ${q.fraction_of_pdb_ci95[s].join(" to ")})`).join("; ")
        + ". Fingerprints are KiSSim-lite.";
    }
  })();

  // ---------------------------------------------------------------- Q6
  K.json("/api/summary/induced_fit").then((rows) => {
    const el = document.getElementById("fig-q6");
    const done = (rows || []).filter((r) => r.status === "done");
    if (!done.length) return pending(el, "The induced-fit folds");
    const W = 640, H = 40 + done.length * 22, mg = { t: 20, r: 20, b: 36, l: 170 };
    const x = d3.scaleLinear().domain([0, 1]).range([mg.l, W - mg.r]);
    const y = d3.scaleBand().domain(done.map((r) => `${r.domain_id} · ${r.ccd}`)).range([mg.t, H - mg.b]).padding(0.4);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
    svg.append("g").attr("transform", `translate(0,${H - mg.b})`).call(d3.axisBottom(x).ticks(5).tickFormat(d3.format(".0%")));
    svg.append("g").attr("transform", `translate(${mg.l},0)`).call(d3.axisLeft(y)).selectAll("text").attr("font-family", "JetBrains Mono, monospace").attr("font-size", 10);
    svg.append("text").attr("x", (mg.l + W) / 2).attr("y", H - 4).attr("text-anchor", "middle").attr("font-size", 12).text("Fraction of samples in DFG-out (open: apo, filled: with ligand)");
    done.forEach((r) => {
      const cy = y(`${r.domain_id} · ${r.ccd}`) + y.bandwidth() / 2;
      svg.append("line").attr("x1", x(r.apo_frac_dfgout)).attr("x2", x(r.holo_frac_dfgout)).attr("y1", cy).attr("y2", cy).attr("stroke", C.rule).attr("stroke-width", 3);
      svg.append("circle").attr("cx", x(r.apo_frac_dfgout)).attr("cy", cy).attr("r", 5).attr("fill", "#fff").attr("stroke", C.ink).attr("stroke-width", 1.5);
      svg.append("path").attr("d", `M${x(r.holo_frac_dfgout)},${cy - 6}l5.5,10h-11z`).attr("fill", C.esm).attr("stroke", C.ink)
        .append("title").text(`${r.pdb_id}: apo ${Math.round(r.apo_frac_dfgout * 100)}%, holo ${Math.round(r.holo_frac_dfgout * 100)}% DFG-out; ligand RMSD ${r.lig_rmsd ?? "n/a"} Å`);
    });
  }).catch(() => pending(document.getElementById("fig-q6"), "The induced-fit folds"));
})();
