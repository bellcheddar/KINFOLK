"""Phase 3 (part 2): agreement classes, ensemble summaries, Prior Collapse Index,
pseudokinase motif-repair test (spec 5.6 to 5.9)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import DATA, klifs
from .annotate import core_fit
from .config import log
from .mapping import map_chain
from .structure import kabsch, read_structure

ANN = DATA / "annotation"
STATES = ("BLAminus", "BLAplus", "ABAminus", "BLBminus", "BLBplus", "BLBtrans", "BABtrans", "BBAminus", "other")
RNG = np.random.default_rng(20260927)


def collapse(state: str) -> str:
    return state if state in STATES else "other"


def jsd(p: np.ndarray, q: np.ndarray) -> float:
    p, q = p / p.sum(), q / q.sum()
    m = (p + q) / 2

    def kl(a, b):
        nz = a > 0
        return float(np.sum(a[nz] * np.log2(a[nz] / b[nz])))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def dist(states: list[str]) -> np.ndarray:
    v = np.array([sum(collapse(s) == k for s in states) for k in STATES], dtype=float)
    return v / v.sum() if v.sum() else v


# ------------------------------------------------------------------ agreement
def classify(model_state: str, rmsd: float | None, pdb_states: set[str], model_states: set[str]) -> str:
    """Spec 5.6. `model_states` is the set the model ever samples (ensemble for ESMFold2)."""
    if not pdb_states:
        return "none"
    if model_state not in pdb_states or rmsd is None or not np.isfinite(rmsd) or rmsd > 3.0:
        return "disagree"
    if rmsd >= 1.5:
        return "partial"
    if not pdb_states <= model_states:
        return "partial"
    return "agree"


def entropy(p: np.ndarray) -> float:
    nz = p[p > 0]
    return float(-(nz * np.log2(nz)).sum())


def domain_table(s: pd.DataFrame, cat: pd.DataFrame) -> pd.DataFrame:
    ok = s[s.status == "done"]
    rows = []
    for _, d in cat.iterrows():
        pdb = ok[(ok.source == "pdb") & (ok.domain_id == d.domain_id)]
        af = ok[(ok.source == "afdb") & (ok.domain_id == d.domain_id)]
        es = ok[(ok.source == "esm2") & (ok.domain_id == d.domain_id) & (ok.get("tag", "apo").fillna("apo") == "apo")]
        pdb_states = set(pdb.state) - {"Unassigned"}
        r = {"domain_id": d.domain_id, "n_pdb_chains_q": len(pdb), "n_pdb_states": len(pdb_states),
             "pdb_states": sorted(pdb_states), "pdb_pip": "full" if len(pdb_states) >= 2 else
             "partial" if len(pdb_states) == 1 else "hollow"}
        if len(af):
            a = af.iloc[0]
            r.update(afdb_state=a.state, afdb_rmsd=a.pocket_rmsd_nearest,
                     afdb_agreement=classify(a.state, a.pocket_rmsd_nearest, pdb_states, {a.state}))
        else:
            r.update(afdb_state=None, afdb_rmsd=None, afdb_agreement="none" if not pdb_states else "missing")
        if len(es):
            best = es.sort_values("mean_plddt", ascending=False).iloc[0]
            p = dist(list(es.state))
            ens_states = set(es.state)
            r.update(esm_n=len(es), esm_best=best.struct_id, esm_best_state=best.state,
                     esm_best_rmsd=best.pocket_rmsd_nearest, esm_entropy=round(entropy(p), 3),
                     esm_state_fracs=json.dumps({k: round(float(v), 3) for k, v in zip(STATES, p) if v > 0}),
                     esm_agreement=classify(best.state, best.pocket_rmsd_nearest, pdb_states, ens_states),
                     ensemble_rescue=bool(pdb_states and ((ens_states & pdb_states) - {r.get("afdb_state")})),
                     esm_afdb_rmsd=None)
        else:
            r.update(esm_n=0, esm_agreement="pending")
        rows.append(r)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ PCI
def pci(s: pd.DataFrame, cat: pd.DataFrame, dom: pd.DataFrame) -> dict:
    ok = s[s.status == "done"]
    groups = dict(zip(cat.domain_id, cat.group))
    have = dom[(dom.n_pdb_states > 0) & dom.afdb_state.notna() & (dom.esm_n > 0)].domain_id.tolist()
    if len(have) < 5:
        return {"status": "pending", "n_kinases": len(have)}
    per = {}
    for d in have:
        sub = ok[ok.domain_id == d]
        es = sub[(sub.source == "esm2") & (sub.get("tag", "apo").fillna("apo") == "apo")]
        best = dom.set_index("domain_id").loc[d, "esm_best_state"]
        per[d] = {"pdb": dist([x for x in sub[sub.source == "pdb"].state if x != "Unassigned"]),
                  "afdb": dist(list(sub[sub.source == "afdb"].state)),
                  "esm_single": dist([best]),
                  "esm_ensemble": dist(list(es.state))}
    comps = ("afdb", "esm_single", "esm_ensemble")

    def pooled(keys):
        return {src: np.mean([per[k][src] for k in keys], axis=0) for src in ("pdb", *comps)}

    def stat(keys):
        P = pooled(keys)
        return {c: jsd(P["pdb"], P[c]) for c in comps}

    def boot(keys, n=1000):
        keys = list(keys)
        vals = {c: [] for c in comps}
        for _ in range(n):
            smp = RNG.choice(keys, size=len(keys), replace=True)
            for c, v in stat(smp).items():
                vals[c].append(v)
        return {c: [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]
                for c, v in vals.items()}

    out = {"n_kinases": len(have), "states": list(STATES), "kinome": {}, "by_group": {}}
    point, ci = stat(have), boot(have)
    P = pooled(have)
    out["kinome"] = {c: {"jsd": round(point[c], 4), "ci95": ci[c]} for c in comps}
    out["distributions"] = {k: [round(float(x), 4) for x in v] for k, v in P.items()}
    for g in sorted({groups[d] for d in have}):
        keys = [d for d in have if groups[d] == g]
        if len(keys) < 5:
            continue
        pt, c = stat(keys), boot(keys, 1000)
        out["by_group"][g] = {"n": len(keys), **{k: {"jsd": round(pt[k], 4), "ci95": c[k]} for k in comps},
                              "distributions": {k: [round(float(x), 4) for x in v] for k, v in pooled(keys).items()}}
    return out


# ------------------------------------------------------------------ motif repair
CHECK_POS = {"K": klifs.POS_B3_LYS, "HRD": klifs.POS_HRD_ASP, "DFG": klifs.POS_DFG_ASP}


def _ca_cb(chain, resnum):
    for res in chain:
        if res.seqid.num == resnum and not res.seqid.icode.strip():
            ca = res.find_atom("CA", "*")
            cb = res.find_atom("CB", "*") or ca
            if ca is not None:
                return np.array(ca.pos.tolist()), np.array(cb.pos.tolist())
    return None, None


def motif_repair(s: pd.DataFrame, cat: pd.DataFrame, dom: pd.DataFrame) -> pd.DataFrame:
    from .annotate import pocket_array

    ok = s[s.status == "done"]
    act = cat[~cat.is_pseudokinase]
    act = act[act.domain_id.isin(dom[dom.pdb_states.map(lambda x: "BLAminus" in x)].domain_id)]
    rows = []
    for _, p in cat[cat.is_pseudokinase].iterrows():
        pk = p.klifs_pocket_seq
        ident = act.klifs_pocket_seq.map(lambda q: sum(a == b and a not in "_-" for a, b in zip(pk, q)))
        if ident.empty:
            continue
        hom = act.loc[ident.idxmax()]
        reps = ok[(ok.domain_id == hom.domain_id) & (ok.source == "pdb") & (ok.state == "BLAminus")]
        if reps.empty:
            continue
        ref_row = reps.sort_values("resolution").iloc[0]
        ref_st = read_structure(str(DATA.parent / ref_row.file_path))
        ref_ch = next((x for x in ref_st[0] if x.name == ref_row.chain), ref_st[0][0])
        chains = pd.read_parquet(DATA / "pdb" / "chains.parquet")
        pocket_ref = chains[chains.structure_ID == int(ref_row.struct_id.split(":")[1])].pocket.iloc[0]
        ref_map, _ = map_chain(ref_ch, pocket_ref)
        ref_pk = pocket_array(ref_ch, ref_map)
        degraded = [k for k, pos in CHECK_POS.items() if k in (p.pseudo_evidence or "") or
                    (k == "K" and "K→" in (p.pseudo_evidence or ""))]
        uni_map = klifs.map_pocket_to_sequence(pk, p.full_sequence)
        models = ok[(ok.domain_id == p.domain_id) & ok.source.isin(["afdb", "esm2"]) &
                    (ok.get("tag", "apo").fillna("apo") == "apo")]
        best_esm = dom.set_index("domain_id").get("esm_best", pd.Series(dtype=object)).get(p.domain_id)
        pdb_best = ok[(ok.domain_id == p.domain_id) & (ok.source == "pdb")].sort_values("resolution").head(1)
        for _, m in pd.concat([models[(models.source == "afdb") | (models.struct_id == best_esm)], pdb_best]).iterrows():
            st = read_structure(str(DATA.parent / m.file_path))
            ch = next((x for x in st[0] if len(x) > 20), st[0][0])
            if m.source == "pdb":
                kmap, _ = map_chain(ch, chains[chains.structure_ID == int(m.struct_id.split(":")[1])].pocket.iloc[0])
            elif m.source == "esm2":
                cs = int(m.construct_start)
                kmap = {q: u - cs + 1 for q, u in uni_map.items()}
            else:
                kmap = uni_map
            mob = pocket_array(ch, kmap)
            if core_fit(mob, ref_pk) is None:
                continue
            idx = [q - 1 for q in klifs.CORE_POSITIONS if not np.isnan(mob[q - 1, 0]) and not np.isnan(ref_pk[q - 1, 0])]
            R, t = kabsch(mob[idx], ref_pk[idx])
            for k in degraded:
                pos = CHECK_POS[k]
                if pos not in kmap or pos not in ref_map:
                    continue
                ca, cb = _ca_cb(ch, kmap[pos])
                rca, rcb = _ca_cb(ref_ch, ref_map[pos])
                if ca is None or rca is None:
                    continue
                err = (np.linalg.norm(ca @ R.T + t - rca) + np.linalg.norm(cb @ R.T + t - rcb)) / 2
                call = "restored" if err < 1.5 else "partial" if err < 3.0 else "faithful"
                rows.append({"domain_id": p.domain_id, "struct_id": m.struct_id, "source": m.source, "motif": k,
                             "evidence": p.pseudo_evidence, "homologue": hom.domain_id,
                             "homologue_ref": ref_row.struct_id, "error_A": round(float(err), 2), "call": call})
    return pd.DataFrame(rows)


def run() -> dict:
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    s = pd.read_parquet(ANN / "structures.parquet")
    if "tag" not in s:
        s["tag"] = None
    dom = domain_table(s, cat)
    dom.to_parquet(ANN / "domains.parquet", index=False)
    p = pci(s, cat, dom)
    (ANN / "pci.json").write_text(json.dumps(p, indent=1))
    mr = motif_repair(s, cat, dom)
    mr.to_parquet(ANN / "motif_repair.parquet", index=False)
    log(f"summarise: agreement AFDB {dict(dom.afdb_agreement.value_counts())}; "
        f"ESMFold2 {dict(dom.esm_agreement.value_counts())}; PCI {p.get('kinome', p.get('status'))}; "
        f"motif repair rows {len(mr)}")
    return {"pci": p}
