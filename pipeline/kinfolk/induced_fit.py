"""Phase 7: induced-fit test (spec section 8, Q6).

select()  picks up to 30 ESMFold2-blind, DFG-out, type II complexes with a CCD
          ligand (one per kinase first, best resolution) and queues apo and
          holo folds (8 samples each) in pass Q6 of the fold queue.
run()     reads the annotated samples: DFG-out fraction apo vs holo, Kincore
          label shift, ligand heavy-atom RMSD to the crystal after pocket-core
          superposition (atoms matched by CCD atom name).
"""
from __future__ import annotations

import json

import gemmi
import numpy as np
import pandas as pd

from . import DATA
from .annotate import pocket_array
from .config import decision, log
from .foldqueue import EXTRA_JOBS, construct
from .mapping import kincore_reference, load_mapped
from .structure import kabsch, read_structure

IF = DATA / "induced_fit"
N_MAX = 30
SAMPLES = 8
NON_LIGANDS = {"HOH", "SO4", "GOL", "EDO", "PEG", "CL", "NA", "MG", "ZN", "DMS", "ACT", "PO4", "MN", "CA", "K"}


def select() -> pd.DataFrame:
    IF.mkdir(parents=True, exist_ok=True)
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet").set_index("domain_id")
    ch = pd.read_parquet(DATA / "pdb" / "chains.parquet")
    s = pd.read_parquet(DATA / "annotation" / "structures.parquet")
    lab = s[s.source == "pdb"][["struct_id", "kincore_spatial"]]
    ch["struct_id"] = "pdb:" + ch.structure_ID.astype(str)
    ch = ch.merge(lab, on="struct_id", how="left")
    ref = kincore_reference()
    ref["key"] = ref.pdb + ref.chain
    t2 = set(ref[ref.ligand_types.str.contains("Type2", na=False)].key)
    c = ch[ch.passes_quality & ch.blind_esm & (ch.klifs_ligand_type == "II")
           & ((ch.kincore_spatial == "DFGout") | (ch.DFG == "out"))
           & ~ch.ligand.astype(str).isin(["0", "", "None"])].copy()
    c["kincore_type2"] = (c.pdb.str.upper() + c.chain).isin(t2)
    c["res"] = pd.to_numeric(c.resolution_rcsb.fillna(c.resolution), errors="coerce")
    c = c.sort_values(["kincore_type2", "res"], ascending=[False, True])
    first = c.drop_duplicates("domain_id")
    pick = pd.concat([first, c.drop(first.index)]).head(N_MAX)
    log(f"induced fit: {len(c)} candidate blind DFG-out type II chains over {c.domain_id.nunique()} kinases; "
        f"selected {len(pick)}")
    # ABL1 / imatinib anchor
    abl = ch[(ch.domain_id == "ABL1") & (ch.ligand.astype(str) == "STI") & ch.passes_quality]
    post = abl[abl.blind_esm]
    anchor = {"post_cutoff_entries": sorted(post.pdb.str.upper().unique().tolist()),
              "all_entries": sorted(abl.pdb.str.upper().unique().tolist())}
    anchor["status"] = "blind anchor" if len(post) else "memorisation possible"
    if not len(post) and len(abl):
        best = abl.assign(res=pd.to_numeric(abl.resolution_rcsb.fillna(abl.resolution), errors="coerce")) \
            .sort_values("res").head(1)
        best = best.assign(memorisation_possible=True)
        pick = pd.concat([pick.assign(memorisation_possible=False), best])
    else:
        pick = pick.assign(memorisation_possible=False)
    (IF / "abl1_imatinib.json").write_text(json.dumps(anchor, indent=1))
    pick.to_parquet(IF / "selection.parquet", index=False)

    existing = {(j["domain_id"], j["tag"]) for j in
                (json.loads(l) for l in (EXTRA_JOBS.read_text().splitlines() if EXTRA_JOBS.exists() else []) if l.strip())}
    with open(EXTRA_JOBS, "a") as fh:
        for dom in sorted(set(pick.domain_id)):
            if (dom, "q6apo") in existing:
                continue
            seq, cs = construct(cat.loc[dom])
            fh.write(json.dumps({"domain_id": dom, "tag": "q6apo", "pass": "Q6", "sequence": seq,
                                 "construct_start": cs, "n_samples": SAMPLES, "ligands": []}) + "\n")
        for _, r in pick.iterrows():
            tag = f"holo_{r.pdb.upper()}{r.chain}_{r.ligand}"
            if (r.domain_id, tag) in existing:
                continue
            seq, cs = construct(cat.loc[r.domain_id])
            fh.write(json.dumps({"domain_id": r.domain_id, "tag": tag, "pass": "Q6", "sequence": seq,
                                 "construct_start": cs, "n_samples": SAMPLES, "ligands": [str(r.ligand)]}) + "\n")
    decision("Q6 fold order", "Pass Q6 (apo + holo, 8 samples each) runs after the minimum viable set and T2, "
             "before the extra T3 samples", "Strict tier order (T4 last)",
             "Q6 is a headline question; extra T3 samples only refine ensembles", "Move pass Q6 in foldqueue.build_jobs")
    return pick


def _ligand_atoms(chain_or_model, code: str) -> dict[str, np.ndarray]:
    out = {}
    for ch in chain_or_model:
        for res in ch:
            if res.name == code:
                for a in res:
                    if a.element.name != "H":
                        out.setdefault(a.name, np.array(a.pos.tolist()))
                if out:
                    return out
    return out


def _pred_ligand_atoms(model) -> dict[str, np.ndarray]:
    for ch in model:
        for res in ch:
            if gemmi.find_tabulated_residue(res.name) is None or not gemmi.find_tabulated_residue(res.name).is_amino_acid():
                if res.name not in NON_LIGANDS and len(res) > 5:
                    return {a.name: np.array(a.pos.tolist()) for a in res if a.element.name != "H"}
    return {}


def ligand_rmsd(sample_row, crystal_row, dom) -> float | None:
    ch_p, km_p = load_mapped(sample_row, dom)
    ch_x, km_x = load_mapped(crystal_row, dom)
    P, X = pocket_array(ch_p, km_p), pocket_array(ch_x, km_x)
    from .annotate import CORE_IDX
    ok = CORE_IDX[~np.isnan(P[CORE_IDX, 0]) & ~np.isnan(X[CORE_IDX, 0])]
    if len(ok) < 15:
        return None
    R, t = kabsch(P[ok], X[ok])
    st_p = read_structure(str(DATA.parent / sample_row["file_path"]))
    st_x = gemmi.read_structure(str(DATA.parent / crystal_row["file_path"]))
    lp = _pred_ligand_atoms(st_p[0])
    lx = _ligand_atoms(st_x[0], str(crystal_row["ligand"]))
    common = sorted(set(lp) & set(lx))
    if len(common) < 5:
        return None
    a = np.array([lp[n] for n in common]) @ R.T + t
    b = np.array([lx[n] for n in common])
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def run() -> pd.DataFrame:
    sel_p = IF / "selection.parquet"
    pick = pd.read_parquet(sel_p) if sel_p.exists() else select()
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet").set_index("domain_id")
    s = pd.read_parquet(DATA / "annotation" / "structures.parquet")
    s = s[(s.source == "esm2") & (s.status == "done")]
    allp = pd.read_parquet(DATA / "annotation" / "structures.parquet")
    rows = []
    for _, r in pick.iterrows():
        apo = s[(s.domain_id == r.domain_id) & (s.tag == "q6apo")]
        tag = f"holo_{r.pdb.upper()}{r.chain}_{r.ligand}"
        holo = s[(s.domain_id == r.domain_id) & (s.tag == tag)]
        if apo.empty or holo.empty:
            rows.append({"domain_id": r.domain_id, "ccd": r.ligand, "pdb_id": r.pdb.upper(), "status": "pending"})
            continue
        cryst = allp[allp.struct_id == f"pdb:{r.structure_ID}"]
        rms = []
        if len(cryst):
            for _, h in holo.iterrows():
                try:
                    v = ligand_rmsd(h.to_dict(), cryst.iloc[0].to_dict(), cat.loc[r.domain_id].to_dict())
                    if v is not None:
                        rms.append(v)
                except Exception:  # noqa: BLE001
                    continue
        rows.append({"domain_id": r.domain_id, "ccd": r.ligand, "pdb_id": r.pdb.upper(), "chain": r.chain,
                     "release_date": r.release_date, "memorisation_possible": bool(r.get("memorisation_possible")),
                     "apo_frac_dfgout": round(float((apo.kincore_spatial == "DFGout").mean()), 3),
                     "holo_frac_dfgout": round(float((holo.kincore_spatial == "DFGout").mean()), 3),
                     "apo_states": json.dumps(apo.state.value_counts().to_dict()),
                     "holo_states": json.dumps(holo.state.value_counts().to_dict()),
                     "lig_rmsd": round(float(np.median(rms)), 2) if rms else None,
                     "lig_rmsd_min": round(float(np.min(rms)), 2) if rms else None,
                     "n_apo": len(apo), "n_holo": len(holo), "status": "done"})
    out = pd.DataFrame(rows)
    out.to_parquet(IF / "induced_fit.parquet", index=False)
    done = out[out.status == "done"]
    if len(done):
        log(f"induced fit: {len(done)} complexes; mean DFG-out fraction apo {done.apo_frac_dfgout.mean():.2f} "
            f"vs holo {done.holo_frac_dfgout.mean():.2f}; median ligand RMSD {done.lig_rmsd.median()}")
    return out
