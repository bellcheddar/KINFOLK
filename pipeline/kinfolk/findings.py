"""Phase 10: headline findings for Q1 to Q6 -> data/findings.json and FINDINGS.md."""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from . import DATA, ROOT
from .annotate import core_fit, dev_stats, pocket_array
from .config import log, tuning
from .mapping import load_mapped

ANN = DATA / "annotation"


def _j(p, default=None):
    return json.loads(p.read_text()) if p.exists() else default


def _q(v):
    v = pd.Series(v).dropna()
    if not len(v):
        return None
    return {"median": round(float(v.median()), 2), "q1": round(float(v.quantile(.25)), 2),
            "q3": round(float(v.quantile(.75)), 2), "n": int(len(v))}


# ------------------------------------------------------------------ Q3 blind comparison
def _blind_job(args):
    pdb_row, preds, dom = args
    out = []
    try:
        ch, km = load_mapped(pdb_row, dom)
        X = pocket_array(ch, km)
    except Exception:  # noqa: BLE001
        return out
    for pr in preds:
        try:
            c2, k2 = load_mapped(pr, dom)
            fit = core_fit(pocket_array(c2, k2), X)
            if fit is not None:
                out.append({"domain_id": dom["domain_id"], "source": pr["source"], "pdb": pdb_row["struct_id"],
                            "rmsd": float(dev_stats(fit, X)[0])})
        except Exception:  # noqa: BLE001
            continue
    return out


def blind_comparison(s: pd.DataFrame, cat: pd.DataFrame, dom: pd.DataFrame) -> dict:
    """Predicted pocket vs the best-resolution chain of each category, per domain (85-position Ca RMSD)."""
    ok = s[s.status == "done"]
    pdb = ok[ok.source == "pdb"].copy()
    pdb["res"] = pd.to_numeric(pdb.resolution, errors="coerce").fillna(9.9)
    catd = cat.set_index("domain_id")
    best = dom.set_index("domain_id").esm_best.dropna().to_dict() if "esm_best" in dom else {}
    cats = {
        "Training era (before 2018-04-30)": pdb[~pdb.blind_af2.fillna(False).astype(bool)],
        "Blind for ESMFold2 (after 2021-09-30)": pdb[pdb.blind_esm.fillna(False).astype(bool)],
        "Blind and sequence-novel": pdb[pdb.novel_esm.fillna(False).astype(bool)],
    }
    res = {}
    for label, sub in cats.items():
        jobs = []
        for d, g in sub.sort_values("res").groupby("domain_id"):
            preds = []
            af = ok[(ok.struct_id == f"afdb:{d}")]
            if len(af):
                preds.append(af.iloc[0].to_dict())
            if d in best and (ok.struct_id == best[d]).any():
                preds.append(ok[ok.struct_id == best[d]].iloc[0].to_dict())
            if preds:
                jobs.append((g.iloc[0].to_dict(), preds, catd.loc[d].to_dict() | {"domain_id": d}))
        rows = []
        with ProcessPoolExecutor(int(tuning().get("cpu_workers", 4))) as ex:
            for r in ex.map(_blind_job, jobs, chunksize=8):
                rows += r
        df = pd.DataFrame(rows)
        res[label] = {src: _q(df[df.source == src].rmsd) for src in ("afdb", "esm2")} if len(df) else {}
    return res


# ------------------------------------------------------------------ assemble
def run() -> dict:
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    s = pd.read_parquet(ANN / "structures.parquet")
    dom = pd.read_parquet(ANN / "domains.parquet")
    motif = _j(ANN / "motif_summary.json", {})
    pci = _j(ANN / "pci.json", {})
    q5 = _j(DATA / "selectivity" / "q5.json", {})
    ctrl = _j(DATA / "pockets" / "controls.json", [])
    sites = pd.read_parquet(DATA / "pockets" / "sites.parquet") if (DATA / "pockets" / "sites.parquet").exists() else pd.DataFrame()
    ifit = pd.read_parquet(DATA / "induced_fit" / "induced_fit.parquet") if (DATA / "induced_fit" / "induced_fit.parquet").exists() else pd.DataFrame()
    mr = pd.read_parquet(ANN / "motif_repair.parquet") if (ANN / "motif_repair.parquet").exists() else pd.DataFrame()
    F: dict = {}

    # Q1
    a, e = motif.get("afdb", {}), motif.get("esm2_best", {})
    worst_a = max(a, key=lambda m: a[m]["median"]) if a else None
    F["q1"] = {"afdb": a, "esm2_best": e, "figure": "app /almanac#q1",
               "summary": (f"Against the nearest experimental state, AlphaFold DB motifs deviate by a median of "
                           + ", ".join(f"{m} {a[m]['median']:.2f} Å" for m in ("g.l", "aC", "xDFG", "a.l") if m in a)
                           + (f"; ESMFold2's best sample by " + ", ".join(f"{m} {e[m]['median']:.2f} Å" for m in ("g.l", "aC", "xDFG", "a.l") if m in e) if e else "")
                           + f". The least faithful motif is the {worst_a} region, as expected for the most mobile elements.")
               if a else "Pending."}
    # Q2
    if pci.get("kinome"):
        k = pci["kinome"]
        F["q2"] = {**pci, "summary": (
            f"Across {pci['n_kinases']} kinases, the Prior Collapse Index (Jensen–Shannon divergence from the PDB state "
            f"distribution) is {k['afdb']['jsd']:.3f} for AlphaFold DB (95% CI {k['afdb']['ci95'][0]:.3f} to {k['afdb']['ci95'][1]:.3f}), "
            f"{k['esm_single']['jsd']:.3f} for the best ESMFold2 sample ({k['esm_single']['ci95'][0]:.3f} to {k['esm_single']['ci95'][1]:.3f}) "
            f"and {k['esm_ensemble']['jsd']:.3f} for the ESMFold2 ensemble ({k['esm_ensemble']['ci95'][0]:.3f} to {k['esm_ensemble']['ci95'][1]:.3f}). "
            + ("The ensemble recovers more of the PDB's state diversity than a single model."
               if k['esm_ensemble']['jsd'] < min(k['afdb']['jsd'], k['esm_single']['jsd']) else
               "The ensemble does not recover more of the PDB's state diversity than a single model."))}
        if "ensemble_rescue" in dom:
            F["q2"]["ensemble_rescue_n"] = int(dom.ensemble_rescue.fillna(False).astype(bool).sum())
    else:
        F["q2"] = {"summary": "Pending: needs ESMFold2 samples for kinases with PDB structures."}
    # Q3
    blind = blind_comparison(s, cat, dom)
    dark = dom[dom.domain_id.isin(cat[cat.is_dark].domain_id)]
    calls = mr[mr.source.isin(["afdb", "esm2"])].groupby("source").call.value_counts().unstack(fill_value=0).to_dict("index") if len(mr) else {}
    F["q3"] = {"blind": blind, "dark_esm_afdb_rmsd": _q(dark.get("esm_afdb_rmsd", pd.Series(dtype=float))),
               "motif_repair": calls}
    bn = blind.get("Blind and sequence-novel", {})
    tr = blind.get("Training era (before 2018-04-30)", {})
    parts = []
    if bn.get("esm2") and tr.get("esm2"):
        parts.append(f"ESMFold2's median pocket RMSD is {tr['esm2']['median']:.2f} Å on training-era structures and "
                     f"{bn['esm2']['median']:.2f} Å on blind, sequence-novel ones (n {bn['esm2']['n']})")
    if bn.get("afdb") and tr.get("afdb"):
        parts.append(f"AlphaFold DB's is {tr['afdb']['median']:.2f} and {bn['afdb']['median']:.2f} Å respectively")
    if F["q3"]["dark_esm_afdb_rmsd"]:
        parts.append(f"on the dark kinome the two predictors' pockets differ by a median {F['q3']['dark_esm_afdb_rmsd']['median']:.2f} Å")
    for src, name in (("afdb", "AlphaFold DB"), ("esm2", "ESMFold2")):
        c = calls.get(src)
        if c:
            tot = sum(c.values())
            parts.append(f"{name} restores the degraded catalytic residue to canonical geometry in {c.get('restored', 0)} of {tot} pseudokinase motifs")
    F["q3"]["summary"] = ("; ".join(parts) + ".") if parts else "Pending."
    F["q3"]["summary"] = F["q3"]["summary"][0].upper() + F["q3"]["summary"][1:]
    # Q4
    rec = [c for c in ctrl if c.get("recovered")]
    ev = [c for c in ctrl if c.get("recovered") is not None]
    n_cand = int(sites.candidate.sum()) if len(sites) else 0
    af_rec = [c for c in ctrl if c.get("afdb_recovered")]
    n_kwc = int(sites[sites.candidate].domain_id.nunique()) if len(sites) else 0
    F["q4"] = {"controls": ctrl, "n_candidates": n_cand, "n_kinases_with_candidates": n_kwc,
               "summary": (f"Apo ESMFold2 ensembles recover {len(rec)} of {len(ev)} positive-control allosteric sites "
                           f"({', '.join(c['gene'] for c in rec) or 'none'}); the single AlphaFold DB model recovers {len(af_rec)}. "
                           f"Kinome-wide, {n_cand} candidate cryptic sites pass the filters across "
                           f"{n_kwc} kinases.")
               if ev else "Pending."}
    # Q5
    if q5.get("rho"):
        fr = q5.get("fraction_of_pdb", {})
        F["q5"] = {**q5, "summary": (
            f"Over {q5['n_kinases']} kinases in the Davis panel, KiSSim-lite pocket similarity correlates with inhibitor "
            f"co-inhibition at Spearman ρ {q5['rho'].get('pdb', float('nan')):.2f} for PDB pockets "
            f"(95% CI {q5['rho_ci95']['pdb'][0]:.2f} to {q5['rho_ci95']['pdb'][1]:.2f}). "
            + " ".join(f"{'ESMFold2' if s_ == 'esm2' else 'AlphaFold DB'} pockets reach {v:.2f} of the PDB correlation "
                       f"(95% CI {q5['fraction_of_pdb_ci95'][s_][0]:.2f} to {q5['fraction_of_pdb_ci95'][s_][1]:.2f})."
                       for s_, v in fr.items()))}
    else:
        F["q5"] = {"summary": "Pending."}
    # Q6
    done = ifit[ifit.status == "done"] if len(ifit) else pd.DataFrame()
    if len(done):
        blind = done[~done.memorisation_possible.fillna(False).astype(bool)]
        p = wilcoxon(blind.holo_frac_dfgout, blind.apo_frac_dfgout).pvalue if len(blind) >= 6 and \
            (blind.holo_frac_dfgout != blind.apo_frac_dfgout).any() else None
        F["q6"] = {"n": int(len(blind)), "apo_mean": round(float(blind.apo_frac_dfgout.mean()), 3),
                   "holo_mean": round(float(blind.holo_frac_dfgout.mean()), 3),
                   "wilcoxon_p": None if p is None else float(p),
                   "lig_rmsd": _q(blind.lig_rmsd),
                   "summary": (f"On {len(blind)} post-cutoff type II complexes, folding with the ligand raises the share of "
                               f"DFG-out samples from {blind.apo_frac_dfgout.mean():.0%} (apo) to {blind.holo_frac_dfgout.mean():.0%} (holo)"
                               + (f" (paired Wilcoxon p = {p:.2g})" if p is not None else "")
                               + (f"; the median ligand heavy-atom RMSD to the crystal is {blind.lig_rmsd.median():.1f} Å "
                                  f"over the {int(blind.lig_rmsd.notna().sum())} complexes where ligand atoms could be "
                                  f"matched by name" if blind.lig_rmsd.notna().any() else "") + ".")}
    else:
        F["q6"] = {"summary": "Pending."}
    (DATA / "findings.json").write_text(json.dumps(F, indent=1, default=str))
    write_md(F)
    log("findings: written data/findings.json and FINDINGS.md")
    return F


def write_md(F: dict) -> None:
    titles = {"q1": "Q1: Motif fidelity", "q2": "Q2: Conformational prior (Prior Collapse Index)",
              "q3": "Q3: Novel sequences and structures", "q4": "Q4: Cryptic pockets",
              "q5": "Q5: Selectivity", "q6": "Q6: Induced fit"}
    figs = {"q1": "/almanac#q1", "q2": "/almanac#q2 (PCI slider)", "q3": "/almanac#q3, /dark, /black-sheep",
            "q4": "/rooms, /almanac#q4", "q5": "/resemblance, /almanac#q5", "q6": "/almanac#q6"}
    lines = ["# KINFOLK findings", "",
             "Headline results for the six questions in the build spec. Intervals are 95% bootstrap intervals over kinases "
             "(1,000 resamples) unless stated. Figures are live in the app at the paths given.", ""]
    for k, t in titles.items():
        lines += [f"## {t}", "", F.get(k, {}).get("summary", "Pending."), "", f"Figure: `{figs[k]}`", ""]
    q3 = F.get("q3", {}).get("blind", {})
    if q3:
        lines += ["### Q3 detail: pocket Cα RMSD to the best-resolution chain in each set (median, IQR, n)", "",
                  "| Set | AlphaFold DB | ESMFold2 best sample |", "|---|---|---|"]
        for lab, v in q3.items():
            f = lambda x: f"{x['median']:.2f} ({x['q1']:.2f} to {x['q3']:.2f}), n {x['n']}" if x else "n/a"  # noqa: E731
            lines.append(f"| {lab} | {f(v.get('afdb'))} | {f(v.get('esm2'))} |")
        lines.append("")
    lines += ["## Caveats", "",
              "- Pocket fingerprints are KiSSim-lite, a reimplementation following the KiSSim design; kissim itself could not be installed.",
              "- ESMFold2 ran with 3 recycling loops and 100 sampling steps in float32 on Apple MPS (see DECISIONS.md).",
              "- ABL1 with imatinib has no post-cutoff entry, so it is reported separately as \"memorisation possible\".",
              "- PLK5, TAF1_C and TAF1_N could not be placed in KLIFS numbering and are excluded from pocket metrics.",
              "- Neighbour joining needs a complete distance matrix, so the pocket trees also drop domains whose "
              "fingerprints share fewer than 50 KLIFS rows with the rest: EIF2AK4_C from the PDB tree, and PLK5, "
              "STK40, TAF1_C and TAF1_N from the ESMFold2 tree (see data/phylo/pocket_tree_info.json).", ""]
    (ROOT / "FINDINGS.md").write_text("\n".join(lines))
