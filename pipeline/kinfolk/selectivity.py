"""Phase 6: KiSSim-lite pocket fingerprints and Davis 2011 selectivity validation (spec 7.2, Q5).

KiSSim-lite (kissim itself cannot be installed; see DECISIONS.md): for each of
the 85 KLIFS positions, 6 physicochemical features, side-chain exposure
(half-sphere exposure ratio) and Ca distances to the four KiSSim subpocket
centres (hinge region, DFG region, front pocket, pocket centre). Features are
scaled to [0, 1]; distance = mean absolute difference over shared
(position, feature) pairs; similarity = 1 - distance.
"""
from __future__ import annotations

import json
import re
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import DATA
from .config import log, tuning
from .mapping import load_mapped

SEL = DATA / "selectivity"
SIZE = {**dict.fromkeys("GASC", 1), **dict.fromkeys("PTVDNILM", 2), **dict.fromkeys("EQKHFRYW", 3)}
HBD = {"R": 3, "K": 1, "N": 1, "Q": 1, "H": 1, "S": 1, "T": 1, "W": 1, "Y": 1}
HBA = {"D": 2, "E": 2, "N": 1, "Q": 1, "H": 1, "S": 1, "T": 1, "Y": 1}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
AROM = set("FWYH")
ALIPH = set("AILMPV")
CENTRES = {"hinge": (16, 47, 80), "dfg": (19, 24, 81), "front": (6, 48, 75)}
N_FEAT = 11
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
         "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
         "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M", "TPO": "T", "SEP": "S", "PTR": "Y"}


# ------------------------------------------------------------------ fingerprint
def _hse_ratio(ca_all: np.ndarray, ca: np.ndarray, cb: np.ndarray, radius=12.0) -> float:
    v = cb - ca
    if np.linalg.norm(v) < 1e-3:
        return np.nan
    d = ca_all - ca
    near = (np.linalg.norm(d, axis=1) < radius) & (np.linalg.norm(d, axis=1) > 0.1)
    up = (d[near] @ v > 0).sum()
    return up / max(1, near.sum())


def fingerprint(chain, kmap: dict[int, int]) -> np.ndarray:
    res = {r.seqid.num: r for r in chain if not r.seqid.icode.strip()}
    ca_all = np.array([a.pos.tolist() for r in chain for a in r if a.name == "CA"])
    ca = np.full((85, 3), np.nan)
    fp = np.full((85, N_FEAT), np.nan)
    for p, n in kmap.items():
        r = res.get(n)
        if r is None or r.find_atom("CA", "*") is None:
            continue
        ca[p - 1] = r.find_atom("CA", "*").pos.tolist()
    centres = {}
    for k, idx in CENTRES.items():
        pts = ca[[i - 1 for i in idx]]
        pts = pts[~np.isnan(pts[:, 0])]
        centres[k] = pts.mean(0) if len(pts) == len(idx) else None
    good = ~np.isnan(ca[:, 0])
    centres["centre"] = ca[good].mean(0) if good.sum() >= 40 else None
    for p, n in kmap.items():
        r = res.get(n)
        if r is None or np.isnan(ca[p - 1, 0]):
            continue
        aa = THREE.get(r.name, "X")
        cb_atom = r.find_atom("CB", "*")
        if cb_atom is not None:
            cb = np.array(cb_atom.pos.tolist())
        else:  # glycine: virtual CB from backbone
            nn, cc = r.find_atom("N", "*"), r.find_atom("C", "*")
            if nn is None or cc is None:
                cb = ca[p - 1]
            else:
                b = ca[p - 1] - np.array(nn.pos.tolist())
                c = np.array(cc.pos.tolist()) - ca[p - 1]
                cb = -0.58273431 * np.cross(b, c) + 0.56802827 * b - 0.54067466 * c + ca[p - 1]
        f = [SIZE.get(aa, 2) / 3, HBD.get(aa, 0) / 3, HBA.get(aa, 0) / 2, (CHARGE.get(aa, 0) + 1) / 2,
             float(aa in AROM), float(aa in ALIPH), _hse_ratio(ca_all, ca[p - 1], cb)]
        for k in ("hinge", "dfg", "front", "centre"):
            c = centres[k]
            f.append(np.nan if c is None else float(np.clip((np.linalg.norm(ca[p - 1] - c) - 2) / 22, 0, 1)))
        fp[p - 1] = f
    return fp


def fp_distance(a: np.ndarray, b: np.ndarray) -> float:
    ok = ~np.isnan(a) & ~np.isnan(b)
    rows = ok.any(1).sum()
    if rows < 50:
        return np.nan
    return float(np.abs(a[ok] - b[ok]).mean())


def _fp_job(args):
    row, dom = args
    try:
        ch, kmap = load_mapped(row, dom)
        return row["domain_id"], row["source"], fingerprint(ch, kmap)
    except Exception:  # noqa: BLE001
        return row["domain_id"], row["source"], None


def fingerprints() -> dict[str, dict[str, np.ndarray]]:
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet").set_index("domain_id")
    s = pd.read_parquet(DATA / "annotation" / "structures.parquet")
    dom = pd.read_parquet(DATA / "annotation" / "domains.parquet")
    s = s[s.status == "done"]
    jobs = []
    pdb = s[s.source == "pdb"].copy()
    pdb["res"] = pd.to_numeric(pdb.resolution, errors="coerce").fillna(9.9)
    for d, g in pdb.sort_values("res").groupby("domain_id"):
        jobs.append((g.iloc[0].to_dict() | {"source": "pdb"}, cat.loc[d].to_dict()))
    for _, r in s[s.source == "afdb"].iterrows():
        jobs.append((r.to_dict(), cat.loc[r.domain_id].to_dict()))
    best = set(dom.esm_best.dropna()) if "esm_best" in dom else set()
    for _, r in s[s.struct_id.isin(best)].iterrows():
        jobs.append((r.to_dict() | {"source": "esm2"}, cat.loc[r.domain_id].to_dict()))
    out: dict[str, dict[str, np.ndarray]] = {"pdb": {}, "afdb": {}, "esm2": {}}
    with ProcessPoolExecutor(int(tuning().get("cpu_workers", 4))) as ex:
        for d, src, fp in ex.map(_fp_job, jobs, chunksize=8):
            if fp is not None:
                out[src][d] = fp
    SEL.mkdir(parents=True, exist_ok=True)
    for src, v in out.items():
        if v:
            np.savez_compressed(SEL / f"fp_{src}.npz", ids=np.array(list(v)), fps=np.array(list(v.values())))
    log(f"selectivity: KiSSim-lite fingerprints PDB {len(out['pdb'])}, AFDB {len(out['afdb'])}, "
        f"ESMFold2 {len(out['esm2'])}")
    return out


def distance_matrix(fps: dict[str, np.ndarray]) -> tuple[list[str], np.ndarray]:
    names = sorted(fps)
    n = len(names)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = fp_distance(fps[names[i]], fps[names[j]])
    return names, D


# ------------------------------------------------------------------ Davis
JAK_DOM = {"JH1domain-catalytic": "JH1", "JH2domain-pseudokinase": "JH2"}


def davis_hits(cat: pd.DataFrame) -> pd.DataFrame:
    from tdc.multi_pred import DTI

    df = DTI(name="DAVIS", path=str(SEL / "tdc")).get_data()
    ids = set(cat.domain_id)
    genes = cat.groupby("gene").domain_id.apply(list).to_dict()

    def to_domain(t: str) -> str | None:
        m = re.match(r"^([A-Za-z0-9]+)(?:\((.*)\))?$", t.replace("-phosphorylated", "").replace("-nonphosphorylated", ""))
        if not m:
            return None
        gene, var = m.group(1), m.group(2)
        if "phosphorylated" in t or (var and re.search(r"[A-Z]\d+[A-Z]|del|ins|ITD", var)):
            return None  # wild type only
        if var in JAK_DOM and f"{gene}_{JAK_DOM[var]}" in ids:
            return f"{gene}_{JAK_DOM[var]}"
        if var and "Kin.Dom.1" in var:
            return f"{gene}_N" if f"{gene}_N" in ids else None
        if var and "Kin.Dom.2" in var:
            return f"{gene}_C" if f"{gene}_C" in ids else None
        dl = genes.get(gene, [])
        return dl[0] if len(dl) == 1 else None

    df["domain_id"] = df.Target_ID.map(to_domain)
    unmapped = sorted(set(df[df.domain_id.isna()].Target_ID))
    (SEL / "davis_unmapped.txt").write_text("\n".join(unmapped) + "\n")
    df = df.dropna(subset=["domain_id"])
    df["hit"] = df.Y <= 1000.0  # Kd <= 1 uM (TDC reports Kd in nM)
    M = df.pivot_table(index="domain_id", columns="Drug_ID", values="hit", aggfunc="max").fillna(False).astype(bool)
    log(f"selectivity: Davis mapped {len(M)} wild-type domains x {M.shape[1]} inhibitors "
        f"({len(unmapped)} targets unmapped or variant)")
    return M


def coinhibition(M: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    X = M.values.astype(float)
    inter = X @ X.T
    cnt = X.sum(1)
    union = cnt[:, None] + cnt[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        J = np.where(union > 0, inter / union, np.nan)
    return list(M.index), J


def correlate(fp_D: dict[str, tuple[list[str], np.ndarray]], M: pd.DataFrame, n_boot=1000) -> dict:
    ids_c, J = coinhibition(M)
    common = set(ids_c)
    for names, _ in fp_D.values():
        common &= set(names)
    common = sorted(common)
    ci = {d: i for i, d in enumerate(ids_c)}

    fidx = {src: {d: i for i, d in enumerate(names)} for src, (names, _) in fp_D.items()}

    def rho(keys, src):
        D = fp_D[src][1]
        a, b = np.triu_indices(len(keys), 1)
        f = np.array([fidx[src][k] for k in keys])
        c = np.array([ci[k] for k in keys])
        sim = 1 - D[f[a], f[b]]
        co = J[c[a], c[b]]
        ok = ~np.isnan(sim) & ~np.isnan(co)
        return spearmanr(sim[ok], co[ok]).statistic if ok.sum() > 10 else np.nan

    rng = np.random.default_rng(11)
    point = {src: rho(common, src) for src in fp_D}
    boots = {src: [] for src in fp_D}
    ratios = {src: [] for src in fp_D if src != "pdb"}
    for _ in range(n_boot):
        smp = sorted(set(rng.choice(common, size=len(common), replace=True)))
        vals = {src: rho(smp, src) for src in fp_D}
        for src, v in vals.items():
            boots[src].append(v)
            if src != "pdb" and vals.get("pdb"):
                ratios[src].append(v / vals["pdb"])
    q = lambda v: [round(float(np.nanpercentile(v, 2.5)), 3), round(float(np.nanpercentile(v, 97.5)), 3)]  # noqa: E731
    out = {"n_kinases": len(common), "n_pairs": len(common) * (len(common) - 1) // 2,
           "rho": {s: round(float(v), 3) for s, v in point.items()}, "rho_ci95": {s: q(v) for s, v in boots.items()},
           "fraction_of_pdb": {s: round(float(point[s] / point["pdb"]), 3) for s in ratios if point.get("pdb")},
           "fraction_of_pdb_ci95": {s: q(v) for s, v in ratios.items()},
           "note": "Bootstrap over kinases (1,000 resamples; duplicates collapsed). Fingerprints are KiSSim-lite."}
    return out


def run() -> dict:
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    fps = fingerprints()
    fp_D = {src: distance_matrix(v) for src, v in fps.items() if len(v) >= 10}
    for src, (names, D) in fp_D.items():
        np.savez_compressed(SEL / f"dist_{src}.npz", ids=np.array(names), D=D)
    M = davis_hits(cat)
    res = correlate(fp_D, M)
    (SEL / "q5.json").write_text(json.dumps(res, indent=1))
    ids_c, J = coinhibition(M)
    ci = {d: i for i, d in enumerate(ids_c)}
    rows = []
    for src, (names, D) in fp_D.items():
        for i, a in enumerate(names):
            order = np.argsort(np.nan_to_num(D[i], nan=9.0))
            for j in order[1:51]:
                b = names[j]
                rows.append({"domain_a": a, "domain_b": b, "source": src,
                             "kissim_sim": None if np.isnan(D[i, j]) else round(1 - float(D[i, j]), 4),
                             "davis_coinhib": None if a not in ci or b not in ci or np.isnan(J[ci[a], ci[b]])
                             else round(float(J[ci[a], ci[b]]), 4)})
    pd.DataFrame(rows).to_parquet(SEL / "similarity.parquet", index=False)
    log(f"selectivity: Q5 Spearman {res['rho']} (CI {res['rho_ci95']}); fraction of PDB {res['fraction_of_pdb']}")
    # Pocket trees for Phase 4 (PDB and ESMFold2 sources).
    from . import phylo
    phylo.pocket_trees({{"pdb": "pdb", "esm2": "esm"}[k]: v for k, v in fp_D.items() if k in ("pdb", "esm2")})
    return res
