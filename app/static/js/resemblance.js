/* Plate III: tanglegram, off-target finder, Q5 chart. */
(function () {
  "use strict";
  const K = window.KINFOLK;
  const C = K.colours;
  const GROUP_COLOURS = { AGC: "#3F6B52", CAMK: "#4F6C82", CK1: "#8A7A3A", CMGC: "#7C6A92", STE: "#C98F2B",
    TK: "#A9482A", TKL: "#9A5C7A", RGC: "#52625A", Atypical: "#6E8F5B", Other: "#9AA39A" };
  let kin = new Map();

  // ------------------------------------------------------------ tanglegram
  function prune(node, keep) {
    if (!node.children) return keep.has(node.name) ? { ...node } : null;
    const kids = node.children.map((c) => prune(c, keep)).filter(Boolean);
    if (!kids.length) return null;
    if (kids.length === 1) return { ...kids[0], length: (kids[0].length || 0) + (node.length || 0) };
    return { ...node, children: kids };
  }
  function leafNames(n, out = []) { if (!n.children) out.push(n.name); else n.children.forEach((c) => leafNames(c, out)); return out; }

  async function tangle() {
    const src = K.state.get("tg", "pdb");
    document.getElementById("tg-src").value = src;
    const el = document.getElementById("tangle");
    el.innerHTML = "";
    let seq, poc;
    const have = await K.json("/api/trees").catch(() => []);
    if (!have.includes(`tree_pocket_${src}`)) {
      el.innerHTML = `<p class="muted" style="padding:16px">The ${src === "pdb" ? "PDB" : "ESMFold2"} pocket tree is not built yet (Phase 6).</p>`;
      return;
    }
    try {
      [seq, poc] = await Promise.all([K.json("/static/data/tree_sequence.json"), K.json(`/static/data/tree_pocket_${src}.json`)]);
    } catch (e) {
      el.innerHTML = `<p class="muted" style="padding:16px">The ${src === "pdb" ? "PDB" : "ESMFold2"} pocket tree is not available yet.</p>`;
      return;
    }
    const shared = new Set(leafNames(poc));
    const seqP = prune(seq, shared);
    const common = new Set(leafNames(seqP));
    const pocP = prune(poc, common);
    const n = common.size;
    const rowH = 6, W = 980, treeW = 300, H = n * rowH + 40;
    const svg = d3.select(el).append("svg").attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
    const layout = (data, flip) => {
      const h = d3.hierarchy(data, (d) => d.children);
      d3.cluster().size([n * rowH, treeW])(h);
      h.each((d) => { d.X = flip ? W - d.y : d.y; d.Y = d.x + 20; });
      return h;
    };
    const L = layout(seqP, false), R = layout(pocP, true);
    const link = (d) => `M${d.source.X},${d.source.Y}V${d.target.Y}H${d.target.X}`;
    [L, R].forEach((h) => svg.append("g").attr("fill", "none").attr("stroke", C.ink).attr("stroke-width", 0.7)
      .selectAll("path").data(h.links()).enter().append("path").attr("d", link));
    const lpos = new Map(L.leaves().map((l, i) => [l.data.name, [l, i]]));
    const rpos = new Map(R.leaves().map((l, i) => [l.data.name, [l, i]]));
    const pairs = [...common].map((name) => {
      const [a, ia] = lpos.get(name), [b, ib] = rpos.get(name);
      return { name, a, b, shift: Math.abs(ia - ib) / n };
    });
    svg.append("g").selectAll("path").data(pairs).enter().append("path")
      .attr("d", (p) => `M${treeW + 4},${p.a.Y}C${W / 2},${p.a.Y} ${W / 2},${p.b.Y} ${W - treeW - 4},${p.b.Y}`)
      .attr("fill", "none").attr("stroke", (p) => GROUP_COLOURS[(kin.get(p.name) || {}).grp] || "#9AA39A")
      .attr("stroke-width", (p) => (p.shift > 0.25 ? 1.8 : 0.6)).attr("opacity", (p) => (p.shift > 0.25 ? 0.95 : 0.35))
      .append("title").text((p) => `${p.name}: moves ${Math.round(p.shift * 100)}% of the ordering`);
    const disc = pairs.filter((p) => p.shift > 0.25).length;
    const info = K.state.get("tg", "pdb") === "pdb" ? "PDB" : "ESMFold2";
    fetch("/api/summary/tanglegram").then((r) => (r.ok ? r.json() : {})).then((t) => {
      const s = t[src] || {};
      document.getElementById("tg-stats").textContent =
        `${n} kinases shared · ${disc} discordant (shift > 25%)` +
        (s.crossings_before !== undefined ? ` · crossings ${s.crossings_before} → ${s.crossings_after} after untangling` : "") + ` · ${info} pockets`;
    });
  }

  // ------------------------------------------------------------ off-target finder
  async function offTarget() {
    const input = document.getElementById("ot-kinase");
    const srcSel = document.getElementById("ot-src");
    const tbody = document.querySelector("#ot-table tbody");
    const k = input.value.trim();
    if (!kin.has(k)) return;
    K.state.set({ k, src: srcSel.value });
    const rows = (await K.json(`/api/similar/${encodeURIComponent(k)}`)).filter((r) => r.source === srcSel.value);
    tbody.innerHTML = rows.length
      ? rows.slice(0, 25).map((r, i) => `<tr><td>${i + 1}</td><td><a href="/k/${encodeURIComponent(r.domain_b)}">${K.escape(r.domain_b)}</a> <span class="muted">${K.escape((kin.get(r.domain_b) || {}).grp || "")}</span></td>
          <td class="num">${K.fmt(r.kissim_sim, 3)}</td><td class="num">${r.davis_coinhib === null ? "n/a" : K.fmt(r.davis_coinhib, 2)}</td></tr>`).join("")
      : `<tr><td colspan="4" class="muted">No ${srcSel.options[srcSel.selectedIndex].text} fingerprint for ${K.escape(k)}.</td></tr>`;
  }

  // ------------------------------------------------------------ Q5 chart
  function q5() {
    const q = JSON.parse(document.getElementById("q5-data").textContent || "{}");
    const el = document.getElementById("q5-chart");
    if (!q || !q.rho) { el.innerHTML = `<p class="muted">Selectivity analysis pending.</p>`; return; }
    const srcs = [["pdb", "PDB", "circle"], ["afdb", "AlphaFold DB", "square"], ["esm2", "ESMFold2", "triangle"]].filter(([s]) => q.rho[s] !== undefined);
    const W = 460, H = 60 + srcs.length * 54, m = { l: 120, r: 30, t: 20, b: 40 };
    const lo = Math.min(0, ...srcs.map(([s]) => (q.rho_ci95[s] || [0])[0])), hi = Math.max(0.1, ...srcs.map(([s]) => (q.rho_ci95[s] || [0, 0])[1]));
    const x = d3.scaleLinear().domain([lo, hi]).nice().range([m.l, W - m.r]);
    const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
    svg.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(5));
    svg.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", m.t).attr("y2", H - m.b).attr("stroke", C.rule);
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").attr("font-size", 12).text("Spearman ρ with Davis co-inhibition");
    srcs.forEach(([s, label, shape], i) => {
      const y = m.t + 24 + i * 54, [a, b] = q.rho_ci95[s];
      svg.append("text").attr("x", m.l - 10).attr("y", y + 4).attr("text-anchor", "end").attr("font-size", 13).text(label);
      svg.append("line").attr("x1", x(a)).attr("x2", x(b)).attr("y1", y).attr("y2", y).attr("stroke", C.ink).attr("stroke-width", 2);
      const cx = x(q.rho[s]);
      const fill = s === "pdb" ? "#fff" : s === "afdb" ? C.afdb : C.esm;
      if (shape === "circle") svg.append("circle").attr("cx", cx).attr("cy", y).attr("r", 7).attr("fill", fill).attr("stroke", C.ink).attr("stroke-width", 2);
      else if (shape === "square") svg.append("rect").attr("x", cx - 7).attr("y", y - 7).attr("width", 14).attr("height", 14).attr("fill", fill).attr("stroke", C.ink);
      else svg.append("path").attr("d", `M${cx},${y - 8}l7,12h-14z`).attr("fill", fill).attr("stroke", C.ink);
      svg.append("text").attr("x", x(b) + 8).attr("y", y + 4).attr("font-family", "JetBrains Mono, monospace").attr("font-size", 11).text(q.rho[s].toFixed(2));
    });
    const f = q.fraction_of_pdb || {};
    document.getElementById("q5-note").textContent =
      `${q.n_kinases} kinases, ${q.n_pairs} pairs. Predicted-structure correlation as a fraction of the PDB correlation: ` +
      Object.entries(f).map(([s, v]) => `${s === "esm2" ? "ESMFold2" : "AlphaFold DB"} ${v.toFixed(2)} (95% CI ${q.fraction_of_pdb_ci95[s].join(" to ")})`).join("; ") +
      ". Fingerprints are KiSSim-lite.";
  }

  document.getElementById("tg-src").addEventListener("change", (e) => { K.state.set({ tg: e.target.value }); tangle(); });
  document.getElementById("ot-kinase").addEventListener("change", offTarget);
  document.getElementById("ot-src").addEventListener("change", offTarget);
  K.kinases().then((rows) => {
    kin = new Map(rows.map((r) => [r.domain_id, r]));
    document.getElementById("ot-list").innerHTML = rows.map((r) => `<option value="${K.escape(r.domain_id)}">`).join("");
    const k = K.state.get("k", null);
    if (k) { document.getElementById("ot-kinase").value = k; document.getElementById("ot-src").value = K.state.get("src", "pdb"); offTarget(); }
    tangle();
    q5();
  });
})();
