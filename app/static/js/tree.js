/* Plate I: Family Tree (D3 phylogram / cladogram / radial with folded groups). */
(function () {
  "use strict";
  const K = window.KINFOLK;
  const C = K.colours;
  const SPACING = 15;          // px between tips in the unzoomed layout
  const DEPTH = 420;           // px from root to deepest tip (top-down)
  const LABEL_ZOOM = 1.6;      // labels render only above this zoom
  const DUR = 600;
  const STATE_COLOURS = {
    BLAminus: "#3F6B52", BLAplus: "#6E8F5B", ABAminus: "#8A7A3A", BLBminus: "#4F6C82", BLBplus: "#7C6A92",
    BLBtrans: "#9A5C7A", BABtrans: "#C98F2B", BBAminus: "#A9482A", other: "#9AA39A", none: "#FFFFFF",
  };

  const svg = d3.select("#tree");
  const wrap = document.getElementById("tree-wrap");
  const live = document.getElementById("tree-live");
  let root, leaves = [], kin = new Map(), expanded = new Set(), selected = null, focusIdx = -1;
  let zoomK = 1, layout = K.state.get("tree", "phylo"), src = K.state.get("src", "seq"),
    tips = K.state.get("tips", "agree"), support = K.state.get("support", "off");

  const gRoot = svg.append("g");
  const gBands = gRoot.append("g").attr("class", "bands");
  const gLinks = gRoot.append("g").attr("fill", "none").attr("stroke", C.ink).attr("stroke-width", 1.1);
  const gPath = gRoot.append("g").attr("fill", "none").attr("stroke", "#A9482A").attr("stroke-width", 2.6);
  const gNodes = gRoot.append("g");
  const gScale = svg.append("g").attr("class", "scale");
  let tooltip = null;

  const zoom = d3.zoom().scaleExtent([0.2, 12]).on("zoom", (e) => {
    gRoot.attr("transform", e.transform);
    const labelsOn = e.transform.k >= LABEL_ZOOM;
    if (labelsOn !== zoomK >= LABEL_ZOOM) gNodes.selectAll("text.tip").attr("display", labelsOn ? null : "none");
    zoomK = e.transform.k;
    drawScale();
  });
  svg.call(zoom).on("dblclick.zoom", null);

  // ------------------------------------------------------------ data
  function treeUrl() {
    const f = { seq: "tree_sequence", pocket_pdb: "tree_pocket_pdb", pocket_esm: "tree_pocket_esm" }[src] || "tree_sequence";
    return `/static/data/${f}.json`;
  }

  function prepare(data) {
    root = d3.hierarchy(data, (d) => d.children);
    root.each((n) => {
      n.id = n.data.name || `n${Math.random().toString(36).slice(2, 9)}`;
      n.len = (n.parent ? n.parent.len : 0) + (n.parent ? Math.max(0, n.data.length || 0) : 0);
    });
    // Group of each leaf, and folded "pure" clades (all leaves one group, >= 3 leaves).
    root.leaves().forEach((l) => (l.grp = (kin.get(l.data.name) || {}).grp || "Other"));
    // Fold maximal clades that are >= 85% one group (and >= 3 leaves); intruders stay listed in the tooltip.
    root.eachAfter((n) => {
      if (n.children) {
        n.tally = new Map();
        n.children.forEach((c) => c.tally.forEach((v, g) => n.tally.set(g, (n.tally.get(g) || 0) + v)));
        n.count = d3.sum(n.children, (c) => c.count);
      } else {
        n.tally = new Map([[n.grp, 1]]);
        n.count = 1;
      }
      const [g, v] = [...n.tally].sort((a, b) => b[1] - a[1])[0];
      n.pure = v / n.count >= 0.85 ? g : "~mixed";
    });
    root.each((n) => {
      n.foldable = !!n.children && !n.pure.startsWith("~") && n.count >= 3 && (!n.parent || n.parent.pure.startsWith("~"));
    });
  }

  function visibleTree() {
    // Copy of the hierarchy where folded pure clades become single "group" leaves.
    const copy = root.copy();
    copy.each((n) => {
      n.src = origIndex.get(n.data);
      if (n.src && n.src.foldable && !expanded.has(n.src.pure) && !expanded.has(cladeKey(n.src))) {
        n.folded = true;
        n.children = null;
      }
    });
    return copy;
  }
  const origIndex = new Map();
  const cladeKey = (n) => `${n.pure}:${n.leaves()[0].data.name}`;

  // ------------------------------------------------------------ layout
  function compute() {
    root.each((n) => origIndex.set(n.data, n));
    const vt = visibleTree();
    vt.each((n) => {
      n.len = (n.parent ? n.parent.len : 0) + (n.parent ? Math.max(0, n.data.length || 0) : 0);
    });
    const lv = vt.leaves();
    const width = Math.max(600, lv.length * SPACING);
    const lens = lv.map((l) => l.len).sort(d3.ascending);
    const maxLen = Math.max(1e-6, d3.quantile(lens, 0.98) || d3.max(lens) || 1);
    vt.each((n) => { n.clipped = n.len > maxLen * 1.02; });
    const maxDepth = d3.max(lv, (l) => l.depth) || 1;
    d3.cluster().size([width, DEPTH]).separation((a, b) => (a.folded || b.folded ? 2.2 : 1))(vt);
    vt.each((n) => {
      const across = n.x;
      const down = layout === "phylo" ? (Math.min(n.len, maxLen * 1.02) / maxLen) * DEPTH : (n.depth / maxDepth) * DEPTH;
      if (layout === "radial") {
        const a = (across / width) * 2 * Math.PI - Math.PI / 2;
        const r = 40 + (Math.min(n.len, maxLen * 1.02) / maxLen) * DEPTH;
        n.px = Math.cos(a) * r;
        n.py = Math.sin(a) * r;
        n.angle = a;
      } else {
        n.px = across;
        n.py = down;
      }
    });
    vt.scale = maxLen;
    vt.width = width;
    return vt;
  }

  function linkPath(d) {
    const s = d.source, t = d.target;
    if (layout === "radial") return `M${s.px},${s.py}L${t.px},${t.py}`;
    return `M${s.px},${s.py}H${t.px}V${t.py}`;
  }

  // ------------------------------------------------------------ marks
  function pipColour(kind, k) {
    if (!k) return "#fff";
    if (kind === "pdb") return k.pdb_pip === "full" ? C.ink : k.pdb_pip === "partial" ? "#8C9A8F" : "#fff";
    const a = kind === "afdb" ? k.afdb_agreement : k.esm_agreement;
    return { agree: C.agree, partial: C.partial, disagree: C.disagree }[a] || "#fff";
  }

  function drawMarks(sel) {
    sel.selectAll("*").remove();
    sel.each(function (n) {
      const g = d3.select(this);
      if (n.folded) {
        const w = Math.min(60, 8 + n.src.count * 1.2);
        g.append("path").attr("d", `M0,0L${-w / 2},${34}H${w / 2}Z`).attr("fill", "#E9EDE2").attr("stroke", C.ink);
        if (n.src.count >= 8 || zoomK >= LABEL_ZOOM)
          g.append("text").attr("class", "grp-label").attr("transform", "translate(0,42) rotate(90)")
            .attr("font-family", "Instrument Serif, serif").attr("font-style", "italic").attr("font-size", 14)
            .attr("fill", C.ink).text(`${n.src.pure} (${n.src.count})`);
        return;
      }
      const k = kin.get(n.data.name);
      const x0 = 0, y0 = 8;
      if (tips === "agree") {
        const pdb = pipColour("pdb", k);
        g.append("circle").attr("cx", x0).attr("cy", y0).attr("r", 4.2).attr("fill", pdb).attr("stroke", C.ink).attr("stroke-width", 1.2)
          .attr("stroke-dasharray", k && k.pdb_pip === "hollow" ? "2 1.5" : null);
        g.append("rect").attr("x", x0 - 4).attr("y", y0 + 7).attr("width", 8).attr("height", 8).attr("fill", pipColour("afdb", k))
          .attr("stroke", C.ink).attr("stroke-width", 0.8);
        g.append("path").attr("d", `M${x0},${y0 + 18}l4.6,8h-9.2z`).attr("fill", pipColour("esm", k)).attr("stroke", C.ink).attr("stroke-width", 0.8);
      } else if (tips === "state") {
        [["pdb_state", "circle"], ["afdb_state", "rect"], ["esm_state", "tri"]].forEach(([f, shape], i) => {
          const col = STATE_COLOURS[(k && k[f]) || "none"] || STATE_COLOURS.other;
          const y = y0 + i * 9;
          if (shape === "circle") g.append("circle").attr("cx", x0).attr("cy", y).attr("r", 4).attr("fill", col).attr("stroke", C.ink);
          else if (shape === "rect") g.append("rect").attr("x", x0 - 4).attr("y", y - 3).attr("width", 8).attr("height", 8).attr("fill", col).attr("stroke", C.ink).attr("stroke-width", 0.8);
          else g.append("path").attr("d", `M${x0},${y}l4.6,8h-9.2z`).attr("fill", col).attr("stroke", C.ink).attr("stroke-width", 0.8);
        });
      } else if (tips === "cryptic") {
        const n2 = (k && k.cryptic_count) || 0;
        g.append("rect").attr("x", -3).attr("y", y0).attr("width", 6).attr("height", 2 + n2 * 6).attr("fill", n2 ? C.esm : C.rule);
      } else if (tips === "dark") {
        g.append("circle").attr("cx", 0).attr("cy", y0 + 2).attr("r", 5).attr("fill", k && k.is_dark ? C.ink : "#fff").attr("stroke", C.ink);
      }
      const labelY = tips === "agree" || tips === "state" ? y0 + 34 : y0 + 30;
      g.append("text").attr("class", "tip").attr("transform", `translate(0,${labelY}) rotate(90)`)
        .attr("font-family", "JetBrains Mono, monospace").attr("font-size", 9.5).attr("fill", C.ink)
        .attr("display", zoomK >= LABEL_ZOOM ? null : "none").text(n.data.name);
    });
  }

  // ------------------------------------------------------------ render
  let vt = null;
  function render(animate = true) {
    const prev = new Map();
    if (vt) vt.each((n) => prev.set(n.data, [n.px, n.py]));
    vt = compute();
    leaves = vt.leaves();
    const t = svg.transition().duration(animate ? DUR : 0);

    drawBands();
    const links = gLinks.selectAll("path").data(vt.links(), (d) => keyOf(d.target));
    links.exit().remove();
    links.enter().append("path").attr("d", linkPath).attr("opacity", 0).merge(links)
      .attr("stroke-dasharray", (d) => (d.target.clipped && layout === "phylo" ? "4 3" : null))
      .transition(t).attr("opacity", 1).attr("d", linkPath);

    const nodes = gNodes.selectAll("g.node").data(vt.descendants().filter((n) => !n.children || n.folded || support === "on"), (d) => keyOf(d));
    nodes.exit().remove();
    const enter = nodes.enter().append("g").attr("class", "node")
      .attr("transform", (d) => { const p = prev.get(d.data) || [d.px, d.py]; return `translate(${p[0]},${p[1]})`; });
    const merged = enter.merge(nodes);
    merged.on("click", (e, d) => onClick(d)).on("mouseenter", (e, d) => hover(e, d)).on("mouseleave", () => hover(null))
      .style("cursor", (d) => (d.children && !d.folded ? "default" : "pointer"));
    merged.transition(t).attr("transform", (d) => {
      const rot = layout === "radial" && !d.children ? ` rotate(${(d.angle * 180) / Math.PI - 90})` : "";
      return `translate(${d.px},${d.py})${rot}`;
    });
    merged.filter((d) => !d.children || d.folded).call(drawMarks);
    merged.filter((d) => d.children && !d.folded).each(function (d) {
      const g = d3.select(this);
      g.selectAll("*").remove();
      const s = d.data.support;
      if (support === "on" && s !== null && s !== undefined && s < 70)
        g.append("text").attr("x", 3).attr("y", -3).attr("font-family", "JetBrains Mono, monospace").attr("font-size", 8)
          .attr("fill", "#52625A").text(Math.round(s));
    });
    drawPath();
    drawScale();
    document.getElementById("btn-fit-group").disabled = !selected;
  }

  function keyOf(n) {
    return n.data.name || n.data.__k || (n.data.__k = Math.random().toString(36).slice(2));
  }

  function drawBands() {
    gBands.selectAll("*").remove();
    if (layout === "radial") return;
    const runs = [];
    leaves.forEach((l) => {
      const grp = l.folded ? l.src.pure : l.data.name && kin.get(l.data.name) ? kin.get(l.data.name).grp : "Other";
      const last = runs[runs.length - 1];
      if (last && last.grp === grp) last.x1 = l.px;
      else runs.push({ grp, x0: l.px, x1: l.px });
    });
    const sel = selected ? (kin.get(selected) || {}).grp : null;
    gBands.selectAll("rect").data(runs.filter((r) => r.x1 - r.x0 > SPACING * 2)).enter().append("rect")
      .attr("x", (r) => r.x0 - SPACING / 2).attr("width", (r) => r.x1 - r.x0 + SPACING).attr("y", -20).attr("height", DEPTH + 220)
      .attr("fill", (r) => (r.grp === sel ? "#F2E3DA" : "#E9EDE2")).attr("opacity", 0.7);
    gBands.selectAll("text").data(runs.filter((r) => r.x1 - r.x0 > SPACING * 6)).enter().append("text")
      .attr("x", (r) => (r.x0 + r.x1) / 2).attr("y", -6).attr("text-anchor", "middle")
      .attr("font-family", "Instrument Serif, serif").attr("font-style", "italic").attr("font-size", 16).attr("fill", "#52625A")
      .text((r) => r.grp);
  }

  function drawPath() {
    gPath.selectAll("*").remove();
    if (!selected || !vt) return;
    const leaf = leaves.find((l) => l.data.name === selected);
    if (!leaf) return;
    const pts = leaf.ancestors();
    const d = pts.slice(0, -1).map((n) => linkPath({ source: n.parent, target: n })).join("");
    gPath.append("path").attr("d", d);
  }

  function drawScale() {
    gScale.selectAll("*").remove();
    if (layout !== "phylo" || !vt) return;
    const t = d3.zoomTransform(svg.node());
    const px = (0.1 / vt.scale) * DEPTH * t.k;
    const h = svg.node().clientHeight;
    gScale.attr("transform", `translate(16,${h - 24})`);
    gScale.append("line").attr("x1", 0).attr("x2", 0).attr("y1", 0).attr("y2", -px).attr("stroke", C.ink).attr("stroke-width", 2);
    gScale.append("text").attr("x", 6).attr("y", -px / 2).attr("font-family", "JetBrains Mono, monospace").attr("font-size", 11)
      .attr("fill", C.ink).text("0.1 substitutions per site");
  }

  // ------------------------------------------------------------ interaction
  function onClick(d) {
    if (d.folded) {
      expanded.add(cladeKey(d.src));
      K.state.set({ open: [...expanded].join(",") });
      render();
      announce(`Opened a ${d.src.pure} clade, ${d.src.count} kinases`);
      return;
    }
    if (!d.children) select(d.data.name, true);
  }

  function select(name, showTip) {
    selected = name;
    K.state.set({ k: name });
    drawBands();
    drawPath();
    document.getElementById("btn-fit-group").disabled = false;
    const leaf = leaves.find((l) => l.data.name === name);
    if (leaf && showTip) card(leaf);
    dock(name);
    const k = kin.get(name) || {};
    announce(`${name}: ${k.grp || ""}, AlphaFold DB ${K.agreementLabel[k.afdb_agreement] || "n/a"}, ESMFold2 ${K.agreementLabel[k.esm_agreement] || "n/a"}`);
  }

  function card(leaf) {
    if (tooltip) tooltip.remove();
    const k = kin.get(leaf.data.name) || {};
    const t = d3.zoomTransform(svg.node());
    const [x, y] = t.apply([leaf.px, leaf.py]);
    tooltip = d3.select(wrap).append("div").attr("class", "tree-tooltip").attr("role", "dialog").attr("aria-label", `${leaf.data.name} summary`);
    tooltip.html(`<h3>${K.escape(leaf.data.name)}</h3>
      <div class="taxon">${K.escape(k.grp || "")} › ${K.escape(k.family || "")}</div>
      <dl>
        <dt>Nearest by sequence</dt><dd class="mono">${K.escape(k.nearest_seq || "n/a")}</dd>
        <dt>Nearest by pocket</dt><dd class="mono">${K.escape(k.nearest_pocket || "n/a")}</dd>
        <dt>PDB states</dt><dd>${k.n_pdb_states ?? 0}</dd>
        <dt>AlphaFold DB</dt><dd>${K.escape(K.agreementLabel[k.afdb_agreement] || "n/a")}</dd>
        <dt>ESMFold2</dt><dd>${K.escape(K.agreementLabel[k.esm_agreement] || "n/a")}</dd>
      </dl>
      <a class="btn primary" href="/k/${encodeURIComponent(leaf.data.name)}">Open specimen</a>
      <button class="btn" type="button" data-close>Close</button>`);
    const w = wrap.clientWidth;
    tooltip.style("left", `${Math.max(10, Math.min(w - 330, x + 12))}px`).style("top", `${Math.max(10, Math.min(wrap.clientHeight - 260, y + 12))}px`);
    tooltip.select("[data-close]").on("click", () => { tooltip.remove(); tooltip = null; svg.node().focus(); });
  }

  let hoverTip = null;
  function hover(e, d) {
    if (hoverTip) { hoverTip.remove(); hoverTip = null; }
    if (!e || !d || zoomK >= LABEL_ZOOM) return;
    const label = d.folded
      ? `${d.src.pure} clade, ${d.src.count} kinases` + ([...d.src.tally].length > 1
        ? ` (with ${[...d.src.tally].filter(([g]) => g !== d.src.pure).map(([g, v]) => `${v} ${g}`).join(", ")})` : "") + ": click to open"
      : d.data.name;
    if (!label) return;
    const [mx, my] = d3.pointer(e, wrap);
    hoverTip = d3.select(wrap).append("div").attr("class", "tree-tooltip").style("pointer-events", "none")
      .style("left", `${mx + 12}px`).style("top", `${my + 12}px`).style("padding", "6px 10px").text(label);
  }

  function announce(msg) { live.textContent = msg; }

  function dock(name) {
    const sec = document.getElementById("docked");
    K.json(`/static/data/kinase/${encodeURIComponent(name)}.json`).then((d) => {
      sec.hidden = false;
      document.getElementById("docked-title").textContent = `${d.gene} · ${d.domain_id}`;
      document.getElementById("docked-open").href = `/k/${encodeURIComponent(name)}`;
      document.getElementById("docked-body").innerHTML = `
        <p class="taxon">${K.escape(d.group)} › ${K.escape(d.family)} › <span class="mono">${K.escape(d.uniprot)}</span> › <span class="mono">${d.dom_start}-${d.dom_end}</span></p>
        <p class="note">${K.escape(d.observer_note || "")}</p>`;
    }).catch(() => { sec.hidden = true; });
  }

  function fit(nodes) {
    const w = svg.node().clientWidth, h = svg.node().clientHeight;
    const xs = nodes.map((n) => n.px), ys = nodes.map((n) => n.py);
    const x0 = d3.min(xs) - 30, x1 = d3.max(xs) + 30, y0 = d3.min(ys) - 40, y1 = d3.max(ys) + (layout === "radial" ? 40 : 160);
    const k = Math.min(12, 0.95 / Math.max((x1 - x0) / w, (y1 - y0) / h));
    svg.transition().duration(DUR).call(zoom.transform, d3.zoomIdentity.translate(w / 2 - (k * (x0 + x1)) / 2, h / 2 - (k * (y0 + y1)) / 2).scale(k));
  }

  function fitGroup() {
    const g = selected && kin.get(selected) ? kin.get(selected).grp : null;
    if (!g) return;
    if (!expanded.has(g)) { expanded.add(g); render(false); }
    fit(leaves.filter((l) => !l.folded && kin.get(l.data.name) && kin.get(l.data.name).grp === g));
  }

  svg.on("keydown", (e) => {
    const focusables = leaves;
    if (!focusables.length) return;
    if (["ArrowRight", "ArrowDown", "ArrowLeft", "ArrowUp"].includes(e.key)) {
      e.preventDefault();
      focusIdx = (focusIdx + (e.key === "ArrowRight" || e.key === "ArrowDown" ? 1 : -1) + focusables.length) % focusables.length;
      const n = focusables[focusIdx];
      gNodes.selectAll("g.node").attr("outline", null).select(".focus-ring").remove();
      gNodes.selectAll("g.node").filter((d) => d === n).append("circle").attr("class", "focus-ring")
        .attr("r", 11).attr("cy", 8).attr("fill", "none").attr("stroke", "#E07A55").attr("stroke-width", 2.5);
      announce(n.folded ? `${n.src.pure} group, ${n.src.count} kinases, folded` : n.data.name);
    } else if (e.key === "Enter" && focusIdx >= 0) {
      e.preventDefault();
      onClick(focusables[focusIdx]);
    } else if (e.key === "+" || e.key === "=") svg.transition().call(zoom.scaleBy, 1.4);
    else if (e.key === "-") svg.transition().call(zoom.scaleBy, 1 / 1.4);
  });

  // ------------------------------------------------------------ controls
  function segment(id, value, onChange) {
    const el = document.getElementById(id);
    const btns = el.querySelectorAll("button");
    const setPressed = (v) => btns.forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.v === v)));
    setPressed(value);
    btns.forEach((b) => b.addEventListener("click", () => { setPressed(b.dataset.v); onChange(b.dataset.v); }));
  }
  segment("ctl-layout", layout, (v) => { layout = v; K.state.set({ tree: v }); render(); setTimeout(() => fit(vt.descendants()), DUR + 20); });
  segment("ctl-support", support, (v) => { support = v; K.state.set({ support: v === "on" ? "on" : null }); render(false); });
  const srcSel = document.getElementById("ctl-src");
  srcSel.value = src;
  srcSel.addEventListener("change", () => {
    if (srcSel.value === "tangle") { location.href = "/resemblance"; return; }
    src = srcSel.value;
    K.state.set({ src });
    load();
  });
  const tipSel = document.getElementById("ctl-tips");
  tipSel.value = tips;
  tipSel.addEventListener("change", () => { tips = tipSel.value; K.state.set({ tips }); render(false); legend(); });
  document.getElementById("btn-fit").addEventListener("click", () => fit(vt.descendants()));
  document.getElementById("btn-fit-group").addEventListener("click", fitGroup);
  document.getElementById("btn-expand").addEventListener("click", (e) => {
    const all = expanded.size === 0 || e.currentTarget.dataset.state !== "open";
    if (all) root.each((n) => n.foldable && expanded.add(n.pure));
    else expanded.clear();
    e.currentTarget.dataset.state = all ? "open" : "folded";
    e.currentTarget.textContent = all ? "Fold groups" : "Open all";
    K.state.set({ open: [...expanded].join(",") || null });
    render();
    setTimeout(() => fit(vt.descendants()), DUR + 20);
  });

  function legend() {
    const el = document.getElementById("tree-legend");
    if (tips === "state")
      el.innerHTML = Object.entries(STATE_COLOURS).filter(([k]) => k !== "none")
        .map(([k, c]) => `<span><svg width="12" height="12"><rect width="12" height="12" fill="${c}"/></svg>${k}</span>`).join("") +
        "<span>ring PDB (most common) · square AlphaFold DB · triangle ESMFold2 best sample</span>";
    else if (tips === "cryptic") el.innerHTML = "<span>Bar length: candidate cryptic sites from the ESMFold2 ensemble</span>";
    else if (tips === "dark") el.innerHTML = "<span>Filled dot: dark kinome (no PDB chain covering 70% of the domain)</span>";
    else el.innerHTML = "<span>Ring: PDB states (filled 2 or more, grey 1, dashed none) · square AlphaFold DB · triangle ESMFold2 (green agrees, amber partial, rust disagrees, white no reference)</span>";
  }

  function load() {
    K.json(treeUrl()).then((data) => {
      prepare(data);
      render(false);
      const k = K.state.get("k", null);
      if (k && kin.has(k)) {
        const g = kin.get(k).grp;
        expanded.add(g);
        render(false);
        select(k, false);
      }
      fit(vt.descendants());
    });
  }

  K.json("/api/trees").then((have) => {
    [["pocket_pdb", "tree_pocket_pdb"], ["pocket_esm", "tree_pocket_esm"]].forEach(([v, id]) => {
      const opt = srcSel.querySelector(`option[value="${v}"]`);
      if (opt && !have.includes(id)) { opt.disabled = true; opt.textContent += " (pending)"; }
    });
  }).catch(() => {});
  (K.state.get("open", "") || "").split(",").filter(Boolean).forEach((g) => expanded.add(g));
  K.kinases().then((rows) => {
    kin = new Map(rows.map((r) => [r.domain_id, r]));
    legend();
    load();
  });
})();
