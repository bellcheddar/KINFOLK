"""Phase 3: per-structure annotation (spec section 5).

Per domain (one worker process each):
  PDB  every quality chain: KLIFS map, pocket Ca (85 x 3), Kincore labels;
       greedy 1.0 A clustering of pocket Ca -> state representatives
       (highest resolution first).
  AFDB and each ESMFold2 sample: KLIFS map, Kincore, pocket-core
       superposition onto every representative, per-motif deviation (minimum
       over representatives), 85-position ribbon against the nearest
       representative, pLDDT per motif, N-lobe/C-lobe PAE.
Outputs (data/annotation/): structures.parquet, motif_dev.parquet,
ribbon.parquet, pdb_pocket/<domain>.npz. ESMFold2 samples are annotated as the
fold manifest reports them, so this stage can rerun to backfill.
"""
from __future__ import annotations

import gzip
import json
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

from . import DATA, klifs
from .config import Manifest, log, tuning
from .kincore import label as kincore_label
from .mapping import map_chain
from .structure import kabsch, read_structure

ANN = DATA / "annotation"
CORE_IDX = np.array([p - 1 for p in klifs.CORE_POSITIONS])
MOTIF_IDX = {m: np.arange(klifs.REGIONS[m][0] - 1, klifs.REGIONS[m][1]) for m in klifs.MOTIFS}
CLUSTER_RMSD = 1.0
MIN_CORE = 15


def state_of(spatial: str, dihedral: str) -> str:
    """Single categorical state used for distributions and agreement."""
    if dihedral and dihedral != "Unassigned":
        return dihedral
    if spatial and spatial != "Unassigned":
        return f"{spatial}-other"
    return "Unassigned"


# ------------------------------------------------------------------ geometry
def pocket_array(chain, kmap: dict[int, int]) -> np.ndarray:
    arr = np.full((85, 3), np.nan)
    ca = {}
    for res in chain:
        if not res.seqid.icode.strip():
            a = res.find_atom("CA", "*")
            if a is not None:
                ca.setdefault(res.seqid.num, a.pos.tolist())
    for p, n in kmap.items():
        if n in ca:
            arr[p - 1] = ca[n]
    return arr


def core_fit(mob: np.ndarray, ref: np.ndarray) -> np.ndarray | None:
    """Superpose `mob` onto `ref` on the shared pocket-core Ca; return transformed mob."""
    ok = CORE_IDX[~np.isnan(mob[CORE_IDX, 0]) & ~np.isnan(ref[CORE_IDX, 0])]
    if len(ok) < MIN_CORE:
        return None
    R, t = kabsch(mob[ok], ref[ok])
    return mob @ R.T + t


def dev_stats(fit: np.ndarray, ref: np.ndarray) -> tuple[float, np.ndarray]:
    d = np.linalg.norm(fit - ref, axis=1)
    both = ~np.isnan(d)
    return (float(np.sqrt(np.mean(d[both] ** 2))) if both.sum() >= MIN_CORE else np.nan), d


def cluster(pockets: np.ndarray, order: np.ndarray) -> tuple[list[int], np.ndarray]:
    """Greedy leader clustering at CLUSTER_RMSD (all-85 RMSD in the core frame)."""
    reps: list[int] = []
    member = np.full(len(pockets), -1)
    for i in order:
        for ci, r in enumerate(reps):
            fit = core_fit(pockets[i], pockets[r])
            if fit is not None and dev_stats(fit, pockets[r])[0] <= CLUSTER_RMSD:
                member[i] = ci
                break
        else:
            reps.append(int(i))
            member[i] = len(reps) - 1
    return reps, member


# ------------------------------------------------------------------ per domain
def plddt_from_bfactor(chain, kmap) -> dict[int, float]:
    out = {}
    for res in chain:
        a = res.find_atom("CA", "*")
        if a is not None:
            out[res.seqid.num] = a.b_iso
    return out


def lobe_pae(pae: np.ndarray | None, kmap: dict[int, int], first: int) -> float | None:
    """Mean PAE between N-lobe (up to the hinge) and C-lobe residues. Rows indexed resnum-first."""
    if pae is None or pae.ndim != 2 or not kmap:
        return None
    hinge = kmap.get(48) or kmap.get(47) or kmap.get(46)
    lo, hi = min(kmap.values()), max(kmap.values())
    if hinge is None:
        return None
    n = np.arange(lo, hinge + 1) - first
    c = np.arange(hinge + 1, hi + 1) - first
    n, c = n[(n >= 0) & (n < len(pae))], c[(c >= 0) & (c < len(pae))]
    if not len(n) or not len(c):
        return None
    return float((pae[np.ix_(n, c)].mean() + pae[np.ix_(c, n)].mean()) / 2)


def annotate_predicted(struct_id, source, chain, kmap, plddt_by_res, pae, pae_first, reps, rep_pockets,
                       rep_states, extra) -> tuple[dict, list[dict], list[dict]]:
    kc = kincore_label(chain, kmap)
    pk = pocket_array(chain, kmap)
    row = {"struct_id": struct_id, "source": source, **extra,
           "kincore_spatial": kc.spatial, "kincore_dihedral": kc.dihedral,
           "state": state_of(kc.spatial, kc.dihedral), "ac_in": kc.ac_in, "d1": kc.d1, "d2": kc.d2,
           "lys_glu_salt": kc.lys_glu_salt,
           "mean_plddt": float(np.mean([plddt_by_res[n] for n in kmap.values() if n in plddt_by_res]))
           if plddt_by_res else None,
           "pae_lobes": lobe_pae(pae, kmap, pae_first)}
    motifs, ribbon = [], []
    best = (np.inf, None, None)
    motif_min = {m: np.inf for m in klifs.MOTIFS}
    for ri, rp in enumerate(rep_pockets):
        fit = core_fit(pk, rp)
        if fit is None:
            continue
        r, d = dev_stats(fit, rp)
        if np.isfinite(r) and r < best[0]:
            best = (r, ri, d)
        for m, idx in MOTIF_IDX.items():
            dm = d[idx]
            dm = dm[~np.isnan(dm)]
            if len(dm):
                motif_min[m] = min(motif_min[m], float(np.sqrt(np.mean(dm ** 2))))
    row["pocket_rmsd_nearest"] = None if not np.isfinite(best[0]) else round(best[0], 3)
    row["nearest_rep"] = None if best[1] is None else reps[best[1]]
    row["nearest_rep_state"] = None if best[1] is None else rep_states[best[1]]
    inv = {p: n for p, n in kmap.items()}
    for m, idx in MOTIF_IDX.items():
        pl = [plddt_by_res[inv[p + 1]] for p in idx if (p + 1) in inv and inv[p + 1] in plddt_by_res]
        motifs.append({"struct_id": struct_id, "motif": m,
                       "rmsd": None if not np.isfinite(motif_min[m]) else round(motif_min[m], 3),
                       "mean_plddt": round(float(np.mean(pl)), 2) if pl else None})
    for p in range(1, 86):
        n = kmap.get(p)
        dev = None if best[2] is None or np.isnan(best[2][p - 1]) else round(float(best[2][p - 1]), 3)
        ribbon.append({"struct_id": struct_id, "klifs_pos": p, "residue": n, "dev_A": dev,
                       "plddt": round(float(plddt_by_res[n]), 1) if n is not None and n in plddt_by_res else None})
    return row, motifs, ribbon


def process_domain(d: dict, chains: list[dict], esm_jobs: list[dict]) -> dict:
    """Annotate one domain. Returns rows for the three tables plus PDB pocket arrays."""
    rows, motifs, ribbons = [], [], []
    # ---- PDB chains
    pdb_ids, pdb_pockets, pdb_states, pdb_res = [], [], [], []
    for c in chains:
        path = DATA / "pdb" / "raw" / f"{c['structure_ID']}.pdb.gz"
        if not path.exists():
            continue
        try:
            st = read_structure(str(path))
            ch = next((x for x in st[0] if x.name == c["chain"]), st[0][0])
            kmap, ident = map_chain(ch, c["pocket"])
            if len(kmap) < 40:
                continue
            kc = kincore_label(ch, kmap)
            pk = pocket_array(ch, kmap)
        except Exception as e:  # noqa: BLE001
            rows.append({"struct_id": f"pdb:{c['structure_ID']}", "source": "pdb", "domain_id": d["domain_id"],
                         "status": f"failed:{type(e).__name__}"})
            continue
        sid = f"pdb:{c['structure_ID']}"
        state = state_of(kc.spatial, kc.dihedral)
        rows.append({"struct_id": sid, "source": "pdb", "domain_id": d["domain_id"], "pdb_id": c["pdb"].upper(),
                     "chain": c["chain"], "sample_idx": None, "release_date": c.get("release_date"),
                     "resolution": c.get("res"), "blind_af2": c.get("blind_af2"), "blind_esm": c.get("blind_esm"),
                     "novel_af2": c.get("novel_af2"), "novel_esm": c.get("novel_esm"),
                     "kincore_spatial": kc.spatial, "kincore_dihedral": kc.dihedral, "state": state,
                     "ac_in": kc.ac_in, "d1": kc.d1, "d2": kc.d2, "lys_glu_salt": kc.lys_glu_salt,
                     "klifs_dfg": c.get("DFG"), "klifs_ac": c.get("aC_helix"), "ligand": c.get("ligand"),
                     "allosteric_ligand": c.get("allosteric_ligand"), "klifs_ligand_type": c.get("klifs_ligand_type"),
                     "pocket_identity": round(ident, 3), "file_path": str(path.relative_to(DATA.parent)),
                     "status": "done"})
        pdb_ids.append(sid)
        pdb_pockets.append(pk)
        pdb_states.append(state)
        pdb_res.append(float(c.get("res") or 9.9))
    reps: list[str] = []
    rep_pockets: list[np.ndarray] = []
    rep_states: list[str] = []
    if pdb_ids:
        P = np.array(pdb_pockets)
        rep_idx, member = cluster(P, np.argsort(pdb_res))
        reps = [pdb_ids[i] for i in rep_idx]
        rep_pockets = [P[i] for i in rep_idx]
        rep_states = [pdb_states[i] for i in rep_idx]
        for r in rows:
            if r["struct_id"] in pdb_ids:
                i = pdb_ids.index(r["struct_id"])
                r["cluster"] = int(member[i])
                r["is_rep"] = i in rep_idx
        (ANN / "pdb_pocket").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(ANN / "pdb_pocket" / f"{d['domain_id']}.npz", ids=np.array(pdb_ids), pockets=P,
                            states=np.array(pdb_states), reps=np.array(reps))

    uni_map = klifs.map_pocket_to_sequence(d["klifs_pocket_seq"], d["full_sequence"])
    # ---- AFDB
    af = DATA / "afdb" / "raw" / f"{d['uniprot']}.cif.gz"
    if af.exists():
        try:
            st = read_structure(str(af))
            ch = st[0][0]
            pae = None
            pj = DATA / "afdb" / "raw" / f"{d['uniprot']}_pae.json.gz"
            if pj.exists():
                j = json.loads(gzip.decompress(pj.read_bytes()))
                j = j[0] if isinstance(j, list) else j
                pae = np.array(j.get("predicted_aligned_error") or j.get("distance"), dtype=np.float32)
            row, m, rb = annotate_predicted(f"afdb:{d['domain_id']}", "afdb", ch, uni_map,
                                            plddt_from_bfactor(ch, uni_map), pae, 1, reps, rep_pockets, rep_states,
                                            {"domain_id": d["domain_id"], "sample_idx": None, "ptm": None,
                                             "file_path": str(af.relative_to(DATA.parent)), "status": "done"})
            rows.append(row); motifs += m; ribbons += rb
        except Exception as e:  # noqa: BLE001
            rows.append({"struct_id": f"afdb:{d['domain_id']}", "source": "afdb", "domain_id": d["domain_id"],
                         "status": f"failed:{type(e).__name__}: {str(e)[:80]}"})
    # ---- ESMFold2 samples
    for job in esm_jobs:
        cstart = job["construct_start"]
        kmap = {p: u - cstart + 1 for p, u in uni_map.items()}
        for s in job["samples"]:
            idx = s["sample"]
            stem = DATA / "esm2" / "models" / d["domain_id"] / job["tag"] / f"s{idx:02d}"
            sid = f"esm2:{d['domain_id']}:{job['tag']}:{idx}"
            try:
                st = read_structure(f"{stem}.cif.gz")
                ch = next(x for x in st[0] if len(x) > 20)
                npz = np.load(f"{stem}.npz")
                pl = npz["plddt"]
                plddt_by_res = {i + 1: float(v) for i, v in enumerate(pl)}
                pae = npz["pae"].astype(np.float32) if npz["pae"].size else None
                row, m, rb = annotate_predicted(sid, "esm2", ch, kmap, plddt_by_res, pae, 1, reps, rep_pockets,
                                                rep_states, {"domain_id": d["domain_id"], "sample_idx": idx,
                                                             "tag": job["tag"], "ligands": ",".join(job.get("ligands") or []),
                                                             "construct_start": cstart,
                                                             "ptm": s.get("ptm"),
                                                             "file_path": str(stem.relative_to(DATA.parent)) + ".cif.gz",
                                                             "status": "done"})
                rows.append(row); motifs += m; ribbons += rb
            except Exception as e:  # noqa: BLE001
                rows.append({"struct_id": sid, "source": "esm2", "domain_id": d["domain_id"],
                             "status": f"failed:{type(e).__name__}: {str(e)[:80]}"})
    return {"domain_id": d["domain_id"], "rows": rows, "motifs": motifs, "ribbons": ribbons,
            "rep_states": rep_states, "n_reps": len(reps)}


# ------------------------------------------------------------------ driver
def run(domains: list[str] | None = None) -> None:
    ANN.mkdir(parents=True, exist_ok=True)
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    chains = pd.read_parquet(DATA / "pdb" / "chains.parquet")
    chains = chains[chains.passes_quality].copy()
    chains["res"] = pd.to_numeric(chains.resolution_rcsb.fillna(chains.resolution), errors="coerce")
    fold = [r for r in Manifest("fold").records() if r.get("status") == "done"]
    by_dom_fold: dict[str, list] = {}
    for r in fold:
        by_dom_fold.setdefault(r["domain_id"], []).append(r)
    if domains:
        cat = cat[cat.domain_id.isin(domains)]
    parts = {"rows": [], "motifs": [], "ribbons": []}
    workers = int(tuning().get("cpu_workers", 4))
    with ProcessPoolExecutor(workers) as ex:
        futs = {}
        for _, d in cat.iterrows():
            ch = chains[chains.domain_id == d.domain_id].to_dict("records")
            futs[ex.submit(process_domain, d.to_dict(), ch, by_dom_fold.get(d.domain_id, []))] = d.domain_id
        for i, f in enumerate(as_completed(futs), 1):
            try:
                res = f.result()
                for k in parts:
                    parts[k] += res[k]
            except Exception:  # noqa: BLE001
                log(f"annotate: {futs[f]} failed: {traceback.format_exc(limit=1).splitlines()[-1]}")
            if i % 50 == 0:
                log(f"annotate: {i}/{len(futs)} domains")
    s = pd.DataFrame(parts["rows"])
    m = pd.DataFrame(parts["motifs"])
    rb = pd.DataFrame(parts["ribbons"])
    if domains and (ANN / "structures.parquet").exists():
        old = pd.read_parquet(ANN / "structures.parquet")
        keep = ~old.domain_id.isin(domains)
        s = pd.concat([old[keep], s], ignore_index=True)
        oldm, oldr = pd.read_parquet(ANN / "motif_dev.parquet"), pd.read_parquet(ANN / "ribbon.parquet")
        m = pd.concat([oldm[oldm.struct_id.isin(s.struct_id) & ~oldm.struct_id.isin(set(m.struct_id))], m])
        rb = pd.concat([oldr[oldr.struct_id.isin(s.struct_id) & ~oldr.struct_id.isin(set(rb.struct_id))], rb])
    s.to_parquet(ANN / "structures.parquet", index=False)
    m.to_parquet(ANN / "motif_dev.parquet", index=False)
    rb.to_parquet(ANN / "ribbon.parquet", index=False)
    log(f"annotate: {len(s)} structures ({dict(s.source.value_counts())}), "
        f"{int(s.status.astype(str).str.startswith('failed').sum())} failed")
