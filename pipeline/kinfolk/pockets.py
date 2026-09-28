"""Phase 5: pocket detection, site clustering, cryptic candidates and positive controls (spec 7.1, 12).

Residues everywhere are in UniProt numbering (sequence alignment of each chain
to the canonical sequence), so pockets from PDB, AFDB and ESMFold2 compare.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import requests

from . import DATA, ROOT, klifs
from .config import Manifest, decision, log, tuning
from .mapping import chain_sequence

PK = DATA / "pockets"
P2RANK = ROOT / "tools" / "vendor" / "p2rank_2.5.1" / "prank"
PAD = 10
MAX_PDB_PER_DOMAIN = 20
CONTROLS = {"ABL1": "Myristoyl pocket (C-lobe)", "MAP2K1": "Type III pocket adjacent to ATP",
            "PDPK1": "PIF pocket (N-lobe)", "EGFR": "EAI045 allosteric site", "AKT1": "PH-kinase interdomain pocket"}
SUBPOCKET = {}
for _r in ("I", "g.l", "II", "hinge", "linker", "aD", "aE", "VI", "c.l", "VII"):
    SUBPOCKET[_r] = "front"
for _r in ("III", "GK", "VIII"):
    SUBPOCKET[_r] = "gate"
for _r in ("IV", "V"):
    SUBPOCKET[_r] = "back pocket BP-I"
for _r in ("aC", "b.l", "xDFG", "a.l"):
    SUBPOCKET[_r] = "back pocket BP-II"


# ------------------------------------------------------------------ numbering
def uniprot_numbering(chain: gemmi.Chain, full_seq: str) -> dict[int, int]:
    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import SubstitutionMatrix, align_optimal

    seq, nums = chain_sequence(chain)
    if not seq:
        return {}
    clean = lambda s: "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in s)  # noqa: E731
    ali = align_optimal(ProteinSequence(clean(seq)), ProteinSequence(clean(full_seq)),
                        SubstitutionMatrix.std_protein_matrix(), gap_penalty=(-10, -1),
                        terminal_penalty=False, max_number=1)[0]
    return {nums[i]: int(j) + 1 for i, j in ali.trace if i >= 0 and j >= 0}


def write_crop(src: Path, dst: Path, chain_name: str | None, keep: set[int] | None) -> dict[int, int] | None:
    """Write a protein-only PDB of one chain, cropped to residue numbers in `keep` (author numbering)."""
    st = gemmi.read_structure(str(src))
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_hydrogens()
    model = st[0]
    names = [c.name for c in model]
    target = chain_name if chain_name in names else max(model, key=len).name
    for n in names:
        if n != target:
            model.remove_chain(n)
    ch = model[target]
    if keep is not None:
        for i in range(len(ch) - 1, -1, -1):
            if ch[i].seqid.num not in keep:
                del ch[i]
    if len(ch) < 30:
        return None
    ch.name = "A"
    st.setup_entities()
    st.write_pdb(str(dst))
    return {r.seqid.num: r.seqid.num for r in ch}


# ------------------------------------------------------------------ detection
def run_p2rank(files: list[Path], out: Path) -> dict[str, list[dict]]:
    if not P2RANK.exists() or not files:
        return {}
    ds = out / "dataset.ds"
    ds.write_text("\n".join(str(f) for f in files) + "\n")
    subprocess.run([str(P2RANK), "predict", str(ds), "-o", str(out / "p2rank"), "-threads",
                    str(tuning().get("p2rank_threads", 4)), "-visualizations", "0"],
                   check=True, capture_output=True)
    res = {}
    for f in files:
        csv = out / "p2rank" / f"{f.name}_predictions.csv"
        if not csv.exists():
            continue
        df = pd.read_csv(csv, skipinitialspace=True)
        df.columns = [c.strip() for c in df.columns]
        res[f.stem] = [{"rank": int(r["rank"]), "probability": float(r["probability"]),
                        "residues": sorted({int(x.split("_")[1]) for x in str(r["residue_ids"]).split()
                                            if "_" in x and x.split("_")[1].lstrip("-").isdigit()})}
                       for _, r in df.iterrows()]
    return res


def run_fpocket(f: Path, out: Path) -> list[dict]:
    work = out / f.stem
    work.mkdir(exist_ok=True)
    local = work / f.name
    shutil.copy(f, local)
    subprocess.run(["fpocket", "-f", str(local)], check=True, capture_output=True)
    d = work / f"{f.stem}_out"
    info = (d / f"{f.stem}_info.txt").read_text() if (d / f"{f.stem}_info.txt").exists() else ""
    pockets = []
    for block in info.split("Pocket ")[1:]:
        num = int(block.split(":")[0])
        vol = next((float(l.split(":")[1]) for l in block.splitlines() if l.strip().startswith("Volume :")), None)
        atm = d / "pockets" / f"pocket{num}_atm.pdb"
        res = set()
        if atm.exists():
            for line in atm.read_text().splitlines():
                if line.startswith("ATOM"):
                    res.add(int(line[22:26]))
        pockets.append({"rank": num, "volume": vol, "residues": sorted(res)})
    shutil.rmtree(work, ignore_errors=True)
    return pockets


def jaccard(a, b) -> float:
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 0.0


# ------------------------------------------------------------------ per domain
def domain_structures(d: dict, s: pd.DataFrame) -> list[dict]:
    """Structures to scan for one domain: ESMFold2 apo samples, AFDB, PDB reps + apo chains (capped)."""
    out = []
    lo, hi = d["domain_start"] - PAD, d["domain_end"] + PAD
    for _, r in s[(s.domain_id == d["domain_id"]) & (s.status == "done")].iterrows():
        if r.source == "esm2" and (r.get("tag") or "apo") == "apo":
            out.append({"struct_id": r.struct_id, "source": "esm2", "path": r.file_path, "chain": None,
                        "offset": int(r.construct_start) - 1, "crop": None})
        elif r.source == "afdb":
            out.append({"struct_id": r.struct_id, "source": "afdb", "path": r.file_path, "chain": None,
                        "offset": 0, "crop": (lo, hi)})
    pdb = s[(s.domain_id == d["domain_id"]) & (s.source == "pdb") & (s.status == "done")]
    apo = pdb[(pdb.ligand.astype(str).isin(["0", "", "None"])) & (pdb.allosteric_ligand.astype(str).isin(["0", "", "None"]))]
    pick = pd.concat([pdb[pdb.is_rep.fillna(False).astype(bool)], apo]).drop_duplicates("struct_id").head(MAX_PDB_PER_DOMAIN)
    for _, r in pick.iterrows():
        out.append({"struct_id": r.struct_id, "source": "pdb", "path": r.file_path, "chain": r.chain,
                    "offset": None, "crop": None, "apo": r.struct_id in set(apo.struct_id)})
    return out


def scan_domain(d: dict, structs: list[dict]) -> list[dict]:
    """Detect pockets in every structure; return pocket rows with UniProt residue numbers."""
    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        files, meta = [], {}
        for i, x in enumerate(structs):
            src = ROOT / x["path"]
            dst = td / f"s{i}.pdb"
            try:
                st = gemmi.read_structure(str(src))
                ch = next((c for c in st[0] if c.name == x["chain"]), max(st[0], key=len))
                if x["offset"] is not None:
                    to_uni = {r.seqid.num: r.seqid.num + x["offset"] for r in ch}
                else:
                    to_uni = uniprot_numbering(ch, d["full_sequence"])
                keep = None
                if x["crop"]:
                    keep = {n for n, u in to_uni.items() if x["crop"][0] <= u <= x["crop"][1]}
                if write_crop(src, dst, ch.name, keep) is None:
                    continue
                files.append(dst)
                meta[dst.stem] = (x, to_uni)
            except Exception:  # noqa: BLE001
                continue
        p2 = {}
        try:
            p2 = run_p2rank(files, td)
        except Exception:  # noqa: BLE001
            p2 = {}
        for f in files:
            x, to_uni = meta[f.stem]
            try:
                fp = run_fpocket(f, td)
            except Exception:  # noqa: BLE001
                fp = []
            for p in p2.get(f.stem, []):
                res = sorted({to_uni[r] for r in p["residues"] if r in to_uni})
                vol = max((q["volume"] for q in fp if jaccard(q["residues"], p["residues"]) >= 0.3 and q["volume"]),
                          default=None)
                rows.append({"domain_id": d["domain_id"], "struct_id": x["struct_id"], "source": x["source"],
                             "apo": x.get("apo"), "tool": "p2rank", "rank": p["rank"],
                             "probability": p["probability"], "volume": vol, "residues": res})
            if not p2:
                for q in fp:
                    res = sorted({to_uni[r] for r in q["residues"] if r in to_uni})
                    rows.append({"domain_id": d["domain_id"], "struct_id": x["struct_id"], "source": x["source"],
                                 "apo": x.get("apo"), "tool": "fpocket", "rank": q["rank"], "probability": None,
                                 "volume": q["volume"], "residues": res})
    return rows


# ------------------------------------------------------------------ sites
def klifs_positions(d: dict) -> dict[int, int]:
    """UniProt residue -> KLIFS position."""
    m = klifs.map_pocket_to_sequence(d["klifs_pocket_seq"], d["full_sequence"])
    return {u: p for p, u in m.items()}


def subpocket(residues: list[int], kpos: dict[int, int]) -> tuple[str, float]:
    inside = [klifs.region_of(kpos[r]) for r in residues if r in kpos]
    frac = len(inside) / max(1, len(residues))
    if frac < 0.3:
        return "outside KLIFS (allosteric)", frac
    c = Counter(SUBPOCKET[x] for x in inside)
    return c.most_common(1)[0][0], frac


def cluster_sites(d: dict, pockets: pd.DataFrame) -> list[dict]:
    kpos = klifs_positions(d)
    esm = pockets[(pockets.source == "esm2")]
    n_samples = esm.struct_id.nunique()
    sites: list[dict] = []
    for _, p in esm.sort_values("probability", ascending=False).iterrows():
        if len(p.residues) < 4:
            continue
        for s in sites:
            if jaccard(p.residues, s["residues"]) >= 0.5:
                s["members"].append(p)
                break
        else:
            sites.append({"residues": list(p.residues), "members": [p]})
    pdb = pockets[pockets.source == "pdb"]
    pdb_apo = pdb[pdb.apo.fillna(False).astype(bool)]
    af = pockets[pockets.source == "afdb"]
    out = []
    for i, s in enumerate(sites):
        m = pd.DataFrame(s["members"])
        counts = Counter(r for rs in m.residues for r in rs)
        consensus = sorted(r for r, c in counts.items() if c >= max(1, len(m) / 2))
        sp, frac_klifs = subpocket(consensus, kpos)
        in_cleft = sp != "outside KLIFS (allosteric)" and frac_klifs >= 0.5
        known = any(jaccard(consensus, q) >= 0.5 for q in pdb.residues)
        in_apo = any(jaccard(consensus, q) >= 0.5 for q in pdb_apo.residues)
        recovery = m.struct_id.nunique() / max(1, n_samples)
        pmax = float(m.probability.max()) if m.probability.notna().any() else None
        vols = m.volume.dropna()
        out.append({"site_id": f"{d['domain_id']}:S{i + 1}", "domain_id": d["domain_id"], "source": "esm2",
                    "recovery": round(recovery, 3), "p2rank_max": None if pmax is None else round(pmax, 3),
                    "vol_min": float(vols.min()) if len(vols) else None,
                    "vol_max": float(vols.max()) if len(vols) else None,
                    "subpocket": sp, "frac_in_klifs": round(frac_klifs, 2), "atp_cleft": in_cleft,
                    "known_in_pdb": known, "in_pdb_apo": in_apo,
                    "in_afdb": any(jaccard(consensus, q) >= 0.5 for q in af.residues),
                    "candidate": bool(recovery >= 0.25 and not in_apo and (pmax or 0) >= 0.5 and not in_cleft
                                      and not known),
                    "n_samples": n_samples, "residues": consensus,
                    "klifs_positions": [kpos[r] for r in consensus if r in kpos]})
    return out


# ------------------------------------------------------------------ controls
def _rcsb_cif(pdb_id: str) -> Path:
    p = PK / "control_entries" / f"{pdb_id.lower()}.cif.gz"
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id.lower()}.cif.gz", timeout=120)
        r.raise_for_status()
        p.write_bytes(r.content)
    return p


def control_sites(cat: pd.DataFrame) -> list[dict]:
    """Resolve each control site from ligand-bound PDB entries (residues within 4 A of the ligand)."""
    chains = pd.read_parquet(DATA / "pdb" / "chains.parquet")
    out = []
    for gene, site in CONTROLS.items():
        d = cat[cat.gene == gene].iloc[0].to_dict()
        c = chains[(chains.domain_id == d["domain_id"]) & chains.passes_quality &
                   ~chains.allosteric_ligand.astype(str).isin(["0", "", "None"])].copy()
        if c.empty:
            out.append({"gene": gene, "domain_id": d["domain_id"], "site_name": site, "status": "failed:no allosteric chain"})
            continue
        lig = c.allosteric_ligand.astype(str).value_counts().index[0]
        c["res_num"] = pd.to_numeric(c.resolution_rcsb.fillna(c.resolution), errors="coerce")
        best = c[c.allosteric_ligand.astype(str) == lig].sort_values("res_num").iloc[0]
        try:
            st = gemmi.read_structure(str(_rcsb_cif(best.pdb)))
            st.setup_entities()
            model = st[0]
            prot = model[best.chain] if best.chain in [x.name for x in model] else max(model, key=len)
            ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
            lig_atoms = [(res, a) for ch in model for res in ch if res.name == lig for a in res]
            contacts = set()
            for _, a in lig_atoms:
                for mk in ns.find_atoms(a.pos, "\0", radius=4.0):
                    cra = mk.to_cra(model)
                    if cra.chain.name == prot.name and cra.residue.name != lig and cra.residue.het_flag != "H":
                        contacts.add(cra.residue.seqid.num)
            to_uni = uniprot_numbering(prot, d["full_sequence"])
            uni = sorted({to_uni[n] for n in contacts if n in to_uni})
            out.append({"gene": gene, "domain_id": d["domain_id"], "site_name": site, "pdb_id": best.pdb.upper(),
                        "chain": prot.name, "ligand": lig, "residues": uni, "status": "done"})
            log(f"controls: {gene} site from {best.pdb.upper()} {prot.name} ligand {lig}: {len(uni)} residues")
        except Exception as e:  # noqa: BLE001
            out.append({"gene": gene, "domain_id": d["domain_id"], "site_name": site, "status": f"failed:{e}"})
    return out


def score_controls(ctrl: list[dict], pockets: pd.DataFrame) -> list[dict]:
    rows = []
    for c in ctrl:
        if c.get("status") != "done":
            rows.append({**c, "recovered": None, "score": None, "afdb_recovered": None})
            continue
        dom_tags = [c["domain_id"]]
        p = pockets[pockets.domain_id.isin(dom_tags)]
        if c["gene"] == "AKT1":
            p = pockets[(pockets.domain_id == "AKT1") & (pockets.struct_id.str.contains(":phkd:") |
                                                         (pockets.source == "afdb"))]
        esm = p[p.source == "esm2"]
        n = esm.struct_id.nunique()
        best_j = esm.groupby("struct_id").residues.apply(lambda rs: max(jaccard(c["residues"], r) for r in rs)) \
            if n else pd.Series(dtype=float)
        hit = (best_j >= 0.4)
        af = p[p.source == "afdb"]
        af_j = max((jaccard(c["residues"], r) for r in af.residues), default=0.0)
        rows.append({**c, "n_samples": n, "recovered": bool(hit.any()) if n else None,
                     "score": round(float(hit.mean()), 3) if n else None,
                     "best_jaccard": round(float(best_j.max()), 3) if n else None,
                     "afdb_jaccard": round(af_j, 3), "afdb_recovered": af_j >= 0.4})
    return rows


# ------------------------------------------------------------------ driver
def _scan(args):
    d, structs = args
    try:
        return d["domain_id"], scan_domain(d, structs), None
    except Exception as e:  # noqa: BLE001
        return d["domain_id"], [], f"{type(e).__name__}: {e}"


def run() -> dict:
    PK.mkdir(parents=True, exist_ok=True)
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    s = pd.read_parquet(DATA / "annotation" / "structures.parquet")
    for col in ("tag", "construct_start", "is_rep", "ligand", "allosteric_ligand", "chain"):
        if col not in s:
            s[col] = None
    # AKT1 PH + kinase samples are scanned as part of AKT1.
    man = Manifest("pocket_structs")      # one line per scanned structure (resumable, incremental)
    scanned = man.done()
    store = PK / "pockets.parquet"
    old = pd.read_parquet(store) if store.exists() else pd.DataFrame()
    tasks = []
    for _, d in cat.iterrows():
        structs = domain_structures(d.to_dict(), s)
        if d.gene == "AKT1":
            for _, r in s[(s.domain_id == "AKT1") & (s.get("tag") == "phkd") & (s.status == "done")].iterrows():
                structs.append({"struct_id": r.struct_id, "source": "esm2", "path": r.file_path, "chain": None,
                                "offset": int(r.construct_start) - 1, "crop": None})
        todo = [x for x in structs if x["struct_id"] not in scanned]
        if todo:
            tasks.append((d.to_dict(), todo))
    log(f"pockets: {sum(len(t[1]) for t in tasks)} new structures to scan over {len(tasks)} domains")
    rows, n = [], 0
    with ProcessPoolExecutor(max(1, int(tuning().get("cpu_workers", 4)) // 4)) as ex:
        for (d, todo), (dom, r, err) in zip(tasks, ex.map(_scan, tasks)):
            rows += r
            n += 1
            for x in todo:
                man.add(x["struct_id"], status="done" if err is None else f"failed:{err}", domain_id=dom)
            if n % 25 == 0:
                log(f"pockets: scanned {n}/{len(tasks)} domains")
    new = pd.DataFrame(rows)
    pockets = pd.concat([old, new], ignore_index=True) if len(old) else new
    if len(pockets):
        pockets = pockets.drop_duplicates(["struct_id", "tool", "rank"], keep="last")
    pockets.to_parquet(store, index=False)

    sites = []
    for _, d in cat.iterrows():
        p = pockets[pockets.domain_id == d.domain_id]
        if len(p):
            sites += cluster_sites(d.to_dict(), p)
    sites_df = pd.DataFrame(sites)
    sites_df["rank_score"] = sites_df.recovery * sites_df.p2rank_max.fillna(0)
    sites_df.to_parquet(PK / "sites.parquet", index=False)
    ctrl = score_controls(control_sites(cat), pockets)
    (PK / "controls.json").write_text(json.dumps(ctrl, indent=1, default=str))
    rec = [c for c in ctrl if c.get("recovered")]
    log(f"pockets: {len(pockets)} pockets, {len(sites_df)} ESMFold2 sites, "
        f"{int(sites_df.candidate.sum()) if len(sites_df) else 0} candidate cryptic sites; "
        f"controls recovered {len(rec)}/{len(ctrl)}")
    decision("Subpocket assignment", "Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by "
             "majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = "
             ">= 50% inside KLIFS", "KLIFS subpocket centroids (KiFrag)", "Region-based and reproducible from the "
             "KLIFS numbering already computed", "Edit SUBPOCKET in pockets.py")
    return {"sites": len(sites_df), "controls": ctrl}
