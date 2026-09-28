"""Phase 2a/2b: PDB chains via KLIFS (+RCSB metadata) and AlphaFold DB models (spec 4.2, 4.3)."""
from __future__ import annotations

import gzip
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

from . import DATA, klifs
from .config import Manifest, decision, log, tuning

PDB_DIR = DATA / "pdb"
AF_DIR = DATA / "afdb"
CAT = DATA / "catalogue" / "kinases.parquet"
CUTOFF_AF2 = "2018-04-30"
CUTOFF_ESM = "2021-09-30"
BACK_POCKET_II = ("bp_II_out", "bp_II_B", "bp_III", "bp_IV", "bp_V")
BACK_POCKET_I = ("bp_I_A", "bp_I_B", "bp_II_in", "bp_II_A_in", "bp_II_B_in")


@retry(stop=stop_after_attempt(5), wait=wait_exponential(min=2, max=60))
def _get(url: str, **params) -> requests.Response:
    r = requests.get(url, params=params, timeout=120)
    if r.status_code in (404, 400):
        return r
    r.raise_for_status()
    return r


@retry(stop=stop_after_attempt(5), wait=wait_exponential(min=2, max=60))
def _post(url: str, payload: dict) -> dict:
    r = requests.post(url, json=payload, timeout=300)
    r.raise_for_status()
    return r.json()


def ligand_type(s: dict) -> str:
    """KLIFS-subpocket-derived ligand class (I, I1/2, II, III, none)."""
    ortho = s.get("ligand") not in (0, "0", "", None)
    allo = s.get("allosteric_ligand") not in (0, "0", "", None)
    if not ortho:
        return "III" if allo else "none"
    back_ii = any(s.get(k) for k in BACK_POCKET_II)
    back_i = any(s.get(k) for k in BACK_POCKET_I)
    if s.get("DFG") == "out" and back_ii:
        return "II"
    if back_i or back_ii:
        return "I1/2"
    return "I"


# ----------------------------------------------------------------------------- PDB
def klifs_structures(df: pd.DataFrame) -> pd.DataFrame:
    raw = PDB_DIR / "klifs_structures.json"
    if not raw.exists():
        info = json.loads((DATA / "catalogue" / "raw" / "klifs_kinase_information.json").read_text())
        name_to_id = {k["name"]: k["kinase_ID"] for k in info}
        ids = sorted({name_to_id[n] for n in df.klifs_name.dropna() if n in name_to_id})
        out = []
        for i in range(0, len(ids), 40):
            out += klifs.get("structures_list", kinase_ID=",".join(map(str, ids[i:i + 40])))
        raw.write_text(json.dumps(out))
    st = pd.DataFrame(json.loads(raw.read_text()))
    st = st[st.species == "Human"].copy()
    kmap = df.dropna(subset=["klifs_name"]).set_index("klifs_name").domain_id.to_dict()
    st["domain_id"] = st.kinase.map(kmap)
    return st.dropna(subset=["domain_id"])


def rcsb_meta(pdb_ids: list[str]) -> pd.DataFrame:
    cache = PDB_DIR / "rcsb_meta.json"
    meta = json.loads(cache.read_text()) if cache.exists() else {}
    todo = [p for p in pdb_ids if p.upper() not in meta]
    q = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
      rcsb_accession_info{initial_release_date} exptl{method}
      rcsb_entry_info{resolution_combined}}}"""
    for i in range(0, len(todo), 250):
        res = _post("https://data.rcsb.org/graphql", {"query": q, "variables": {"ids": [t.upper() for t in todo[i:i + 250]]}})
        for e in (res.get("data") or {}).get("entries") or []:
            if e is None:
                continue
            meta[e["rcsb_id"]] = {
                "release_date": (e["rcsb_accession_info"] or {}).get("initial_release_date", "")[:10],
                "method": ",".join(x["method"] for x in e.get("exptl") or []),
                "resolution_rcsb": ((e.get("rcsb_entry_info") or {}).get("resolution_combined") or [None])[0],
            }
        cache.write_text(json.dumps(meta))
    return pd.DataFrame.from_dict(meta, orient="index").rename_axis("PDB").reset_index()


def quality(st: pd.DataFrame) -> pd.Series:
    res = pd.to_numeric(st.resolution_rcsb.fillna(st.resolution), errors="coerce")
    xray = st.method.str.contains("X-RAY", na=False) & (res <= 3.0)
    em = st.method.str.contains("ELECTRON MICROSCOPY", na=False) & (res <= 3.5)
    complete = (85 - st.missing_residues.astype(int)) / 85 >= 0.8
    return (xray | em) & complete


def download_klifs_complexes(st: pd.DataFrame) -> None:
    man = Manifest("pdb_download")
    done = man.done()
    raw = PDB_DIR / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    todo = [int(s) for s in st.structure_ID if str(s) not in done]

    def fetch(sid: int):
        try:
            r = _get(f"{klifs.API}/structure_get_pdb_complex", structure_ID=sid)
            if r.status_code != 200 or "ATOM" not in r.text:
                return sid, f"failed:http {r.status_code}"
            with gzip.open(raw / f"{sid}.pdb.gz", "wt") as fh:
                fh.write(r.text)
            return sid, "done"
        except Exception as e:  # noqa: BLE001
            return sid, f"failed:{type(e).__name__}"

    n = 0
    with ThreadPoolExecutor(tuning().get("io_workers", 8)) as ex:
        for sid, status in ex.map(fetch, todo):
            man.add(str(sid), status=status)
            n += 1
            if n % 500 == 0:
                log(f"pdb: downloaded {n}/{len(todo)} KLIFS complexes")


# ------------------------------------------------------------------- novelty (MMseqs2)
def pdb_released_before(date: str) -> set[str]:
    cache = PDB_DIR / f"released_before_{date}.json"
    if cache.exists():
        return set(json.loads(cache.read_text()))
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_accession_info.initial_release_date", "operator": "less_or_equal", "value": f"{date}T23:59:59Z"}},
        "return_type": "entry", "request_options": {"return_all_hits": True}}
    res = _post("https://search.rcsb.org/rcsbsearch/v2/query", q)
    ids = [h["identifier"].upper() for h in res.get("result_set", [])]
    cache.write_text(json.dumps(ids))
    return set(ids)


def seqres() -> dict[str, str]:
    p = PDB_DIR / "pdb_seqres.txt.gz"
    if not p.exists():
        r = requests.get("https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz", timeout=600)
        r.raise_for_status()
        p.write_bytes(r.content)
    out, name = {}, None
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0] if "mol:protein" in line else None
            elif name:
                out[name] = line.strip()
    return out


def novelty(st: pd.DataFrame) -> pd.DataFrame:
    """Max identity of each blind chain to any protein chain released before the cutoff."""
    seqs = seqres()
    work = PDB_DIR / "mmseqs"
    work.mkdir(exist_ok=True)
    for tag, cutoff in (("af2", CUTOFF_AF2), ("esm", CUTOFF_ESM)):
        col = f"maxid_{tag}"
        before = pdb_released_before(cutoff)
        target = work / f"before_{tag}.fa"
        if not target.exists():
            with open(target, "w") as fh:
                for k, v in seqs.items():
                    if k.split("_")[0].upper() in before and len(v) >= 30:
                        fh.write(f">{k}\n{v}\n")
        blind = st[st[f"blind_{tag}"]]
        query = work / f"query_{tag}.fa"
        with open(query, "w") as fh:
            for _, r in blind.iterrows():
                key = f"{r.pdb.lower()}_{r.chain}"
                if key in seqs:
                    fh.write(f">{r.structure_ID}\n{seqs[key]}\n")
        out = work / f"hits_{tag}.m8"
        if not out.exists():
            subprocess.run(["mmseqs", "easy-search", str(query), str(target), str(out), str(work / f"tmp_{tag}"),
                            "--format-output", "query,target,fident,qcov", "-c", "0.5", "--cov-mode", "1",
                            "-e", "1e-3", "--threads", str(tuning().get("cpu_workers", 4))],
                           check=True, capture_output=True)
        hits = pd.read_csv(out, sep="\t", names=["q", "t", "fident", "qcov"]) if out.stat().st_size else \
            pd.DataFrame(columns=["q", "t", "fident", "qcov"])
        best = hits.groupby("q").fident.max().to_dict()
        st[col] = st.structure_ID.map(lambda s: best.get(int(s), best.get(str(s), 0.0)))
        st.loc[~st[f"blind_{tag}"], col] = float("nan")
        st[f"novel_{tag}"] = st[f"blind_{tag}"] & (st[col] < 0.4)
    return st


def acquire_pdb() -> pd.DataFrame:
    df = pd.read_parquet(CAT)
    PDB_DIR.mkdir(exist_ok=True)
    st = klifs_structures(df)
    meta = rcsb_meta(sorted(set(st.pdb.str.upper())))
    st["PDB"] = st.pdb.str.upper()
    st = st.merge(meta, on="PDB", how="left")
    st["passes_quality"] = quality(st)
    st["klifs_ligand_type"] = st.apply(lambda r: ligand_type(r.to_dict()), axis=1)
    st["blind_af2"] = st.release_date > CUTOFF_AF2
    st["blind_esm"] = st.release_date > CUTOFF_ESM
    log(f"pdb: {len(st)} human KLIFS chains, {int(st.passes_quality.sum())} pass the quality filter; "
        f"blind AF2 {int(st.blind_af2.sum())}, blind ESM {int(st.blind_esm.sum())}")
    st = novelty(st)
    log(f"pdb: novel AF2 {int(st.novel_af2.sum())}, novel ESM {int(st.novel_esm.sum())} "
        "(MMseqs2, < 40% identity at >= 50% query coverage)")
    download_klifs_complexes(st[st.passes_quality])
    st.to_parquet(PDB_DIR / "chains.parquet", index=False)
    return st


# ------------------------------------------------------------------------------ AFDB
def acquire_afdb() -> pd.DataFrame:
    df = pd.read_parquet(CAT)
    man = Manifest("afdb")
    done = man.done()
    raw = AF_DIR / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    todo = sorted(set(df.uniprot) - done)

    def fetch(u: str):
        try:
            r = _get(f"https://alphafold.ebi.ac.uk/api/prediction/{u}")
            if r.status_code == 404 or not r.json():
                return u, "absent", {}
            e = r.json()[0]
            cif = _get(e["cifUrl"]).content
            pae = _get(e["paeDocUrl"]).content if e.get("paeDocUrl") else b""
            (raw / f"{u}.cif.gz").write_bytes(gzip.compress(cif))
            if pae:
                (raw / f"{u}_pae.json.gz").write_bytes(gzip.compress(pae))
            return u, "done", {"model_version": e.get("latestVersion"), "entry": e.get("entryId"),
                               "seq_len": len(e.get("uniprotSequence", ""))}
        except Exception as ex:  # noqa: BLE001
            return u, f"failed:{type(ex).__name__}", {}

    with ThreadPoolExecutor(tuning().get("io_workers", 8)) as ex:
        for u, status, extra in ex.map(fetch, todo):
            man.add(u, status=status, **extra)
    recs = {r["id"]: r for r in man.records()}
    df["afdb_status"] = df.uniprot.map(lambda u: recs.get(u, {}).get("status", "failed:missing"))
    covered = df.afdb_status.isin(["done", "absent"]).mean()
    log(f"afdb: {int((df.afdb_status == 'done').sum())} domains with model, "
        f"{int((df.afdb_status == 'absent').sum())} absent; covered {covered:.1%}")
    return df


def assign_t1(st: pd.DataFrame) -> None:
    """T1: domains with at least one blind + novel (ESMFold2) quality chain."""
    df = pd.read_parquet(CAT)
    t1 = set(st[st.passes_quality & st.novel_esm].domain_id)
    mask = df.domain_id.isin(t1) & (df.priority_tier != "T0")
    df.loc[mask, "priority_tier"] = "T1"
    df.to_parquet(CAT, index=False)
    log(f"tiers: {int(mask.sum())} domains promoted to T1 (blind + novel for ESMFold2)")
    if mask.sum() < 3:
        decision("T1 hold-out set is small", f"Only {int(mask.sum())} domains have blind + novel ESMFold2 chains",
                 "Relax novelty to < 50% identity", "Spec defines novel as < 40% identity; kept as specified",
                 "Change the 0.4 threshold in acquire.novelty")
