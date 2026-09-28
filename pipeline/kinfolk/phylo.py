"""Phase 4: sequence tree, pocket trees and tanglegram (spec section 6)."""
from __future__ import annotations

import json
import subprocess
import time
from io import StringIO

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from . import DATA, ROOT
from . import hmm as hmmmod
from .config import decision, log, tuning

PHY = DATA / "phylo"
APP_DATA = ROOT / "app" / "static" / "data"
IQTREE_BUDGET_H = 6


# ------------------------------------------------------------------ sequence tree
def msa() -> tuple[str, int, int]:
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    a2m = hmmmod.hmmalign(dict(zip(cat.domain_id, cat.sequence)))
    cols = {k: "".join(c for c in v if c.isupper() or c == "-") for k, v in a2m.items()}
    L = len(next(iter(cols.values())))
    arr = np.array([list(v) for v in cols.values()])
    keep = (arr == "-").mean(0) <= 0.5
    trimmed = {k: "".join(np.array(list(v))[keep]) for k, v in cols.items()}
    path = PHY / "domains_trimmed.fa"
    path.write_text("".join(f">{k}\n{v}\n" for k, v in trimmed.items()))
    log(f"phylo: MSA {len(trimmed)} x {L} match columns; trimmed {int((~keep).sum())} columns with > 50% gaps")
    return str(path), L, int((~keep).sum())


def run_iqtree(aln: str) -> tuple[str, str]:
    prefix = PHY / "seq_tree"
    tree = prefix.with_suffix(".treefile")
    if tree.exists():
        return str(tree), "cached"
    p = tuning().get("iqtree_threads", "-T AUTO").split()
    t0 = time.time()
    cmd = ["iqtree", "-s", aln, "-m", "LG+G4", "-B", "1000", *p, "--prefix", str(prefix), "-redo", "--quiet"]
    try:
        subprocess.run(cmd, check=True, timeout=IQTREE_BUDGET_H * 3600, capture_output=True)
        method = "iqtree LG+G4 -B 1000"
    except subprocess.TimeoutExpired:
        log("phylo: IQ-TREE ultrafast bootstrap exceeded 6 h; falling back to -fast")
        try:
            subprocess.run(["iqtree", "-s", aln, "-m", "LG+G4", "-fast", *p, "--prefix", str(prefix), "-redo",
                            "--quiet"], check=True, timeout=IQTREE_BUDGET_H * 3600, capture_output=True)
            method = "iqtree LG+G4 -fast"
        except Exception:  # noqa: BLE001
            with open(tree, "w") as fh:
                subprocess.run(["FastTree", "-lg", "-gamma", aln], check=True, stdout=fh)
            method = "FastTree -lg -gamma"
    log(f"phylo: sequence tree by {method} in {(time.time() - t0) / 60:.1f} min")
    return str(tree), method


EPK_GROUPS = ("AGC", "CAMK", "CK1", "CMGC", "STE", "TK", "TKL")


def _purity(tree: Tree, groups: dict[str, str]) -> float:
    """Sum over ePK groups of (members / leaves under their MRCA): 7.0 = all monophyletic."""
    score = 0.0
    for g in EPK_GROUPS:
        leaves = [t for t in tree.get_terminals() if groups.get(t.name) == g]
        if len(leaves) >= 2:
            score += len(leaves) / tree.common_ancestor(leaves).count_terminals()
    return score


def root_tree(tree: Tree, groups: dict[str, str]) -> tuple[Tree, str]:
    """Root on the edge that best separates the ePK groups from the atypicals.

    Every edge is tried as the root; the winner maximises ePK group purity
    (members / leaves under each group's MRCA, summed), ties broken by the
    atypical fraction of the smaller root side.
    """
    import copy

    tree.root_at_midpoint()
    base = _purity(tree, groups)
    best = (base, -1.0, None)
    clades = [c for c in tree.find_clades() if c is not tree.root]
    n = tree.count_terminals()
    for i, c in enumerate(clades):
        t = copy.deepcopy(tree)
        target = list(t.find_clades())[list(tree.find_clades()).index(c)]
        t.root_with_outgroup(target)
        side = target.get_terminals()
        small = side if len(side) <= n / 2 else [x for x in t.get_terminals() if x not in side]
        atyp = sum(groups.get(x.name) not in EPK_GROUPS for x in small) / max(1, len(small))
        s = _purity(t, groups)
        if (s, atyp) > best[:2]:
            best = (s, atyp, i)
    if best[2] is None:
        return tree, f"midpoint (ePK purity {base:.2f} / 7)"
    tree.root_with_outgroup(clades[best[2]])
    return tree, (f"rerooted on the edge maximising ePK group purity ({best[0]:.2f} / 7 vs midpoint {base:.2f}; "
                  f"{best[1]:.0%} non-ePK on the outgroup side)")


def monophyly(tree: Tree, groups: dict[str, str]) -> dict[str, dict]:
    out = {}
    for g in sorted(set(groups.values())):
        leaves = [t for t in tree.get_terminals() if groups.get(t.name) == g]
        if len(leaves) < 2:
            continue
        mrca = tree.common_ancestor(leaves)
        n = mrca.count_terminals()
        intruders = sorted(t.name for t in mrca.get_terminals() if groups.get(t.name) != g)
        out[g] = {"members": len(leaves), "mrca_leaves": n, "monophyletic": n == len(leaves),
                  "intruders": intruders[:25], "n_intruders": len(intruders)}
    return out


# ------------------------------------------------------------------ export
def to_json(clade: Clade, index: dict[str, int]) -> dict:
    node = {"name": clade.name if clade.is_terminal() else None,
            "length": round(float(clade.branch_length or 0.0), 5)}
    if clade.is_terminal():
        node["leaf"] = index.get(clade.name)
    else:
        sup = clade.confidence
        if sup is None and clade.name:
            try:
                sup = float(clade.name)
            except ValueError:
                sup = None
        node["support"] = sup
        node["children"] = [to_json(c, index) for c in clade.clades]
    return node


def export(tree: Tree, name: str, index: dict[str, int]) -> str:
    APP_DATA.mkdir(parents=True, exist_ok=True)
    path = APP_DATA / f"{name}.json"
    path.write_text(json.dumps(to_json(tree.root, index), separators=(",", ":")))
    buf = StringIO()
    Phylo.write(tree, buf, "newick")
    (PHY / f"{name}.nwk").write_text(buf.getvalue())
    return str(path)


def sequence_tree() -> dict:
    PHY.mkdir(parents=True, exist_ok=True)
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    groups = dict(zip(cat.domain_id, cat.group))
    aln, L, trimmed = msa()
    tpath, method = run_iqtree(aln)
    tree = Phylo.read(tpath, "newick")
    tree, choice = root_tree(tree, groups)
    mono = monophyly(tree, groups)
    index = {d: i for i, d in enumerate(cat.domain_id)}
    export(tree, "tree_sequence", index)
    info = {"method": method, "root": choice, "msa_columns": L, "trimmed_columns": trimmed, "monophyly": mono}
    (PHY / "sequence_tree_info.json").write_text(json.dumps(info, indent=1))
    log(f"phylo: rooted {choice}; non-monophyletic groups: "
        + ", ".join(f"{g} ({v['n_intruders']} intruders)" for g, v in mono.items() if not v["monophyletic"]))
    return info


# ------------------------------------------------------------------ NJ + tanglegram
def nj_tree(names: list[str], dist: np.ndarray) -> Tree:
    from skbio import DistanceMatrix
    from skbio.tree import nj

    d = (dist + dist.T) / 2
    np.fill_diagonal(d, 0)
    t = nj(DistanceMatrix(d, ids=names))
    buf = StringIO(str(t))
    tree = Phylo.read(buf, "newick")
    tree.root_at_midpoint()
    return tree


def leaf_order(tree: Tree) -> list[str]:
    return [t.name for t in tree.get_terminals()]


def crossings(order_a: list[str], order_b: list[str]) -> int:
    pos_b = {n: i for i, n in enumerate(order_b)}
    seq = [pos_b[n] for n in order_a if n in pos_b]
    # count inversions (merge sort)
    def sort_count(a):
        if len(a) <= 1:
            return a, 0
        m = len(a) // 2
        l, x = sort_count(a[:m])
        r, y = sort_count(a[m:])
        out, i, j, c = [], 0, 0, x + y
        while i < len(l) and j < len(r):
            if l[i] <= r[j]:
                out.append(l[i]); i += 1
            else:
                out.append(r[j]); j += 1; c += len(l) - i
        out += l[i:] + r[j:]
        return out, c
    return sort_count(seq)[1]


def untangle(fixed: Tree, moving: Tree, passes: int = 3) -> tuple[int, int]:
    """Greedy: at each internal node of `moving`, swap children if it reduces crossings."""
    ref = leaf_order(fixed)
    before = crossings(ref, leaf_order(moving))
    best = before
    for _ in range(passes):
        improved = False
        for clade in moving.get_nonterminals(order="postorder"):
            if len(clade.clades) < 2:
                continue
            clade.clades.reverse()
            c = crossings(ref, leaf_order(moving))
            if c < best:
                best, improved = c, True
            else:
                clade.clades.reverse()
        if not improved:
            break
    return before, best


def pocket_trees(similarity: dict[str, tuple[list[str], np.ndarray]]) -> dict:
    """similarity: source -> (domain_ids, distance matrix). Builds NJ trees and the tanglegram."""
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    index = {d: i for i, d in enumerate(cat.domain_id)}
    info = {}
    for source, (names, dist) in similarity.items():
        t = nj_tree(names, dist)
        keep = set(names)
        seq_sub = Phylo.read(PHY / "tree_sequence.nwk", "newick")
        for leaf in list(seq_sub.get_terminals()):
            if leaf.name not in keep:
                seq_sub.prune(leaf)
        before, after = untangle(seq_sub, t)
        export(t, f"tree_pocket_{source}", index)
        info[source] = {"leaves": len(names), "crossings_before": before, "crossings_after": after}
        log(f"phylo: pocket tree ({source}) {len(names)} leaves; tanglegram crossings {before} -> {after}")
    (PHY / "pocket_tree_info.json").write_text(json.dumps(info, indent=1))
    decision("Tanglegram untangling", "Greedy child swaps on the pocket tree against the fixed sequence tree "
             "(3 passes, accept only swaps that reduce crossings)", "Rotate both trees; exact optimisation",
             "Spec 6.3 asks for a greedy heuristic; the sequence tree stays in its canonical order",
             "Change passes in phylo.untangle")
    return info
