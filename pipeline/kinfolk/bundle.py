"""Phase 8: data bundle for the app (spec 9.1, 10).

app/data/kinfolk.sqlite                  schema of spec 10 (+ helper tables for the UI)
app/static/data/kinase/<domain_id>.json  specimen JSON (<= 150 KB)
app/data/structures/<domain_id>/*.cif.gz domain-cropped mmCIF in the domain's pocket-core frame
"""
from __future__ import annotations

import gzip
import json
import shutil
import sqlite3
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import gemmi
import numpy as np
import pandas as pd

from . import DATA, ROOT, klifs
from .annotate import CORE_IDX, pocket_array
from .config import log, tuning
from .mapping import load_mapped

APP = ROOT / "app"
DB = APP / "data" / "kinfolk.sqlite"
SDIR = APP / "data" / "structures"
JDIR = APP / "static" / "data" / "kinase"
MAX_ESM = 8
PAD = 10
BUDGET_GB = 2.5

SCHEMA = """
CREATE TABLE kinase(domain_id TEXT PRIMARY KEY, gene TEXT, uniprot TEXT, name TEXT, grp TEXT, family TEXT,
  dom_start INT, dom_end INT, is_dark INT, is_pseudo INT, pseudo_evidence TEXT, n_pdb_chains INT,
  n_pdb_states INT, tier TEXT);
CREATE TABLE structure(struct_id TEXT PRIMARY KEY, domain_id TEXT REFERENCES kinase(domain_id),
  source TEXT CHECK(source IN ('pdb','afdb','esm2')), pdb_id TEXT, chain TEXT, sample_idx INT,
  release_date TEXT, resolution REAL, blind_af2 INT, blind_esm INT, novel_af2 INT, novel_esm INT,
  kincore_spatial TEXT, kincore_dihedral TEXT, ac_in INT, d1 REAL, d2 REAL, mean_plddt REAL, ptm REAL,
  pocket_rmsd_nearest REAL, agreement TEXT, file_path TEXT, tag TEXT);
CREATE TABLE motif_dev(struct_id TEXT REFERENCES structure(struct_id), motif TEXT, rmsd REAL, mean_plddt REAL);
CREATE TABLE ribbon(struct_id TEXT REFERENCES structure(struct_id), klifs_pos INT, residue INT, dev_A REAL, plddt REAL);
CREATE TABLE ligand(struct_id TEXT REFERENCES structure(struct_id), ccd TEXT, klifs_type TEXT, in_backpocket INT);
CREATE TABLE site(site_id TEXT PRIMARY KEY, domain_id TEXT REFERENCES kinase(domain_id), source TEXT,
  recovery REAL, p2rank_max REAL, vol_min REAL, vol_max REAL, subpocket TEXT, known_in_pdb INT,
  candidate INT, rank_score REAL, residues_json TEXT);
CREATE TABLE similarity(domain_a TEXT REFERENCES kinase(domain_id), domain_b TEXT REFERENCES kinase(domain_id),
  source TEXT, kissim_sim REAL, davis_coinhib REAL);
CREATE TABLE tree(tree_id TEXT PRIMARY KEY, kind TEXT, newick TEXT, json_path TEXT);
CREATE TABLE control(control_id TEXT PRIMARY KEY, domain_id TEXT REFERENCES kinase(domain_id), site_name TEXT,
  expected_residues_json TEXT, recovered INT, score REAL);
CREATE TABLE induced_fit(domain_id TEXT REFERENCES kinase(domain_id), ccd TEXT, pdb_id TEXT,
  apo_frac_dfgout REAL, holo_frac_dfgout REAL, lig_rmsd REAL);
CREATE TABLE summary(key TEXT PRIMARY KEY, value_json TEXT);
CREATE TABLE kinase_tips(domain_id TEXT PRIMARY KEY, gene TEXT, uniprot TEXT, grp TEXT, family TEXT,
  pdb_pip TEXT, afdb_agreement TEXT, esm_agreement TEXT, pdb_state TEXT, afdb_state TEXT, esm_state TEXT,
  n_pdb_states INT, cryptic_count INT, is_dark INT, is_pseudo INT, nearest_seq TEXT, nearest_pocket TEXT);
CREATE TABLE pseudo_view(domain_id TEXT PRIMARY KEY, gene TEXT, grp TEXT, pseudo_evidence TEXT,
  afdb_call TEXT, esm_call TEXT, pdb_call TEXT, detail_json TEXT, image_struct TEXT);
CREATE TABLE dark_view(domain_id TEXT PRIMARY KEY, gene TEXT, grp TEXT, family TEXT, esm_afdb_rmsd REAL,
  esm_entropy REAL, esm_states TEXT, afdb_plddt REAL, esm_plddt REAL, top_site TEXT, n_candidates INT);
CREATE INDEX ix_structure_domain ON structure(domain_id);
CREATE INDEX ix_structure_domain_source ON structure(domain_id, source);
CREATE INDEX ix_motif_struct ON motif_dev(struct_id);
CREATE INDEX ix_ribbon_struct ON ribbon(struct_id);
CREATE INDEX ix_ligand_struct ON ligand(struct_id);
CREATE INDEX ix_site_domain ON site(domain_id);
CREATE INDEX ix_site_domain_source ON site(domain_id, source);
CREATE INDEX ix_sim_a ON similarity(domain_a);
CREATE INDEX ix_sim_b ON similarity(domain_b);
CREATE INDEX ix_sim_a_source ON similarity(domain_a, source);
CREATE INDEX ix_control_domain ON control(domain_id);
CREATE INDEX ix_if_domain ON induced_fit(domain_id);
"""


def _read(p, **kw) -> pd.DataFrame:
    return pd.read_parquet(p, **kw) if p.exists() else pd.DataFrame()


def _json(p, default):
    return json.loads(p.read_text()) if p.exists() else default


def _clean(o):
    """Recursively replace NaN/inf and numpy scalars so the JSON is strict."""
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def _b(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    return int(bool(v))


def _f(v, d=3):
    try:
        return None if v is None or np.isnan(float(v)) else round(float(v), d)
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------ structures
def frame_for(dom_row: dict, s: pd.DataFrame) -> dict | None:
    """Reference structure for a domain's pocket-core frame: best-resolution PDB rep, else AFDB."""
    pdb = s[(s.source == "pdb") & (s.status == "done")]
    if len(pdb):
        reps = pdb[pdb.is_rep.fillna(False).astype(bool)] if "is_rep" in pdb else pdb
        reps = reps if len(reps) else pdb
        return reps.assign(r=pd.to_numeric(reps.resolution, errors="coerce").fillna(9.9)).sort_values("r").iloc[0].to_dict()
    af = s[(s.source == "afdb") & (s.status == "done")]
    return af.iloc[0].to_dict() if len(af) else None


def esm_keep(es: pd.DataFrame) -> list[str]:
    """Medoid of each Kincore state present (closest to the state's mean d1/d2) + best pLDDT."""
    if es.empty:
        return []
    keep = [es.sort_values("mean_plddt", ascending=False).iloc[0].struct_id]
    for state, g in es.groupby("state"):
        xy = g[["d1", "d2"]].astype(float).fillna(g[["d1", "d2"]].astype(float).mean())
        c = xy.mean().values
        keep.append(g.iloc[int(np.argmin(((xy.values - c) ** 2).sum(1)))].struct_id)
    return list(dict.fromkeys(keep))[:MAX_ESM]


def write_domain_structures(args) -> dict:
    dom, rows, ref = args
    out_dir = SDIR / dom["domain_id"]
    out_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    try:
        ref_ch, ref_map = load_mapped(ref, dom)
        R_pocket = pocket_array(ref_ch, ref_map)
    except Exception:  # noqa: BLE001
        return written
    for r in rows:
        try:
            st = gemmi.read_structure(str(DATA.parent / r["file_path"]))
            st.setup_entities()
            st.remove_hydrogens()
            st.remove_waters()
            ch, kmap = load_mapped(r, dom)
            P = pocket_array(ch, kmap)
            ok = CORE_IDX[~np.isnan(P[CORE_IDX, 0]) & ~np.isnan(R_pocket[CORE_IDX, 0])]
            if len(ok) < 15:
                continue
            from .structure import kabsch
            R, t = kabsch(P[ok], R_pocket[ok])
            tr = gemmi.Transform(gemmi.Mat33(R.tolist()), gemmi.Vec3(*t.tolist()))
            model = st[0]
            keep_name = ch.name
            for c in [c.name for c in model]:
                if c != keep_name and r["source"] != "esm2":
                    model.remove_chain(c)
            # crop protein to the domain +/- PAD (in each file's own numbering)
            if kmap:
                lo, hi = min(kmap.values()) - 25, max(kmap.values()) + 25
                for c in model:
                    for i in range(len(c) - 1, -1, -1):
                        res = c[i]
                        is_poly = gemmi.find_tabulated_residue(res.name) is not None and \
                            gemmi.find_tabulated_residue(res.name).is_amino_acid()
                        if is_poly and not (lo <= res.seqid.num <= hi) and c.name == keep_name:
                            del c[i]
            model.transform_pos_and_adp(tr)
            st.setup_entities()
            doc = st.make_mmcif_document()
            name = r["struct_id"].replace(":", "_") + ".cif.gz"
            with gzip.open(out_dir / name, "wt") as fh:
                fh.write(doc.as_string())
            written[r["struct_id"]] = f"{dom['domain_id']}/{name}"
        except Exception:  # noqa: BLE001
            continue
    return written


# ------------------------------------------------------------------ observer's note
def observer_note(k: dict, d: dict, pdb: pd.DataFrame) -> str:
    af, es = d.get("afdb_state"), d.get("esm_best_state")
    pdb_states = set(d.get("pdb_states") or [])
    parts = []
    if k["is_dark"]:
        rmsd = d.get("esm_afdb_rmsd")
        if af and es and af == es:
            parts.append(f"No experimental structure; both predictors settle this kinase in {af}")
        elif af or es:
            parts.append(f"No experimental structure; AlphaFold DB proposes {af or 'no assignable state'} and "
                         f"ESMFold2 {es or 'no assignable state'}")
        if rmsd is not None:
            parts[-1] += f", with the two pockets {rmsd:.1f} Å apart"
    else:
        if af and es and af == es:
            parts.append(f"Both predictors settle this kinase in {af}")
        else:
            parts.append(f"AlphaFold DB places this kinase in {af or 'no assignable state'}; "
                         f"the best ESMFold2 sample in {es or 'no assignable state'}")
        extra = sorted(pdb_states - {af, es})
        if extra:
            clauses = []
            for sname in extra[:2]:
                lt = set(pdb[pdb.state == sname].klifs_ligand_type.dropna())
                if lt == {"II"}:
                    clauses.append(f"{sname}, seen only with type II ligands")
                elif lt and lt <= {"none"}:
                    clauses.append(f"{sname}, seen only without ligands")
                else:
                    clauses.append(sname)
            parts.append("the PDB also records " + " and ".join(clauses))
        elif pdb_states:
            parts.append("every state in the PDB is among those predicted")
    note = "; ".join(parts) + "."
    if d.get("esm_n", 0) and d.get("esm_state_fracs"):
        fr = json.loads(d["esm_state_fracs"])
        if len(fr) > 1:
            note += f" The ESMFold2 ensemble visits {len(fr)} states ({d.get('esm_entropy', 0):.2f} bits)."
        if d.get("ensemble_rescue"):
            note += " The ensemble reaches a PDB state that the AlphaFold DB model misses."
    if k["is_pseudokinase"]:
        note += f" Pseudokinase ({k['pseudo_evidence']})."
    return note


# ------------------------------------------------------------------ specimen JSON
def field_marks(k: dict) -> dict:
    uni = klifs.map_pocket_to_sequence(k["klifs_pocket_seq"], k["full_sequence"])
    out = {}
    for key, pos in (("b3_lys", klifs.POS_B3_LYS), ("ac_glu", klifs.POS_AC_GLU),
                     ("gatekeeper", klifs.POS_GATEKEEPER), ("dfg_asp", klifs.POS_DFG_ASP)):
        u = uni.get(pos)
        out[key] = {"klifs": pos, "resnum": u, "aa": k["full_sequence"][u - 1] if u else None}
    return out


def specimen_json(k: dict, d: dict, s_dom: pd.DataFrame, fam_pts: list, ribbon: pd.DataFrame, sites: list,
                  sim: pd.DataFrame, nearest_seq: list, files: dict, chains_dom: pd.DataFrame) -> dict:
    pdb = s_dom[(s_dom.source == "pdb") & (s_dom.status == "done")]
    af = s_dom[(s_dom.source == "afdb") & (s_dom.status == "done")]
    es = s_dom[(s_dom.source == "esm2") & (s_dom.status == "done") & (s_dom.tag.fillna("apo") == "apo")]
    rb = {}
    pos_meta = [{"pos": p, "region": klifs.region_of(p), "aa": (k["klifs_pocket_seq"][p - 1]
                 if k["klifs_pocket_seq"][p - 1] not in "_-" else None)} for p in range(1, 86)]
    if len(af):
        r = ribbon[ribbon.struct_id == af.iloc[0].struct_id].sort_values("klifs_pos")
        rb["afdb"] = {"dev": [_f(x, 2) for x in r.dev_A], "residue": [None if pd.isna(x) else int(x) for x in r.residue]}
    best = d.get("esm_best")
    if best:
        r = ribbon[ribbon.struct_id == best].sort_values("klifs_pos")
        rb["esm_best"] = {"dev": [_f(x, 2) for x in r.dev_A], "residue": [None if pd.isna(x) else int(x) for x in r.residue]}
        r = ribbon[ribbon.struct_id.isin(es.struct_id)]
        if len(r):
            m = r.groupby("klifs_pos").dev_A.mean().reindex(range(1, 86))
            rb["esm_mean"] = {"dev": [_f(x, 2) for x in m.values]}
    ligs = chains_dom[~chains_dom.ligand.astype(str).isin(["0", "", "None"])]
    lig_counts = Counter(zip(ligs.ligand.astype(str), ligs.klifs_ligand_type))
    look_pocket = sim[(sim.source == ("pdb" if len(pdb) else "esm2"))].head(6)
    return {
        "domain_id": k["domain_id"], "gene": k["gene"], "name": k["name"], "group": k["group"], "family": k["family"],
        "uniprot": k["uniprot"], "dom_start": int(k["domain_start"]), "dom_end": int(k["domain_end"]),
        "is_dark": bool(k["is_dark"]), "is_pseudo": bool(k["is_pseudokinase"]), "pseudo_evidence": k["pseudo_evidence"],
        "tier": k["priority_tier"], "field_marks": field_marks(k), "observer_note": observer_note(k, d, pdb),
        "klifs_pocket": k["klifs_pocket_seq"],
        "pocket_uniprot": sorted(klifs.map_pocket_to_sequence(k["klifs_pocket_seq"], k["full_sequence"]).values()),
        "agreement": {"pdb_pip": d.get("pdb_pip"), "afdb": d.get("afdb_agreement"), "esm2": d.get("esm_agreement"),
                      "pdb_states": d.get("pdb_states") or [], "afdb_state": d.get("afdb_state"),
                      "esm_best_state": d.get("esm_best_state"), "afdb_rmsd": _f(d.get("afdb_rmsd"), 2),
                      "esm_best_rmsd": _f(d.get("esm_best_rmsd"), 2), "esm_entropy": _f(d.get("esm_entropy"), 3),
                      "esm_state_fracs": json.loads(d["esm_state_fracs"]) if d.get("esm_state_fracs") else {},
                      "ensemble_rescue": bool(d.get("ensemble_rescue")) if d.get("ensemble_rescue") is not None else None},
        "compass": {
            "family": fam_pts,
            "pdb": [{"id": r.struct_id, "pdb": r.pdb_id, "chain": r.chain, "d1": _f(r.d1, 1), "d2": _f(r.d2, 1),
                     "state": r.state, "lig": None if str(r.ligand) in ("0", "None", "") else str(r.ligand),
                     "rep": _b(getattr(r, "is_rep", None)) == 1} for r in pdb.itertuples()][:400],
            "afdb": [{"id": r.struct_id, "d1": _f(r.d1, 1), "d2": _f(r.d2, 1), "state": r.state,
                      "plddt": _f(r.mean_plddt, 1)} for r in af.itertuples()],
            "esm": [{"id": r.struct_id, "i": int(r.sample_idx), "d1": _f(r.d1, 1), "d2": _f(r.d2, 1),
                     "state": r.state, "plddt": _f(r.mean_plddt, 1)} for r in es.itertuples()],
        },
        "ribbon": {"positions": pos_meta, "rows": rb},
        "structures": files,
        "sites": sites[:12],
        "lookalikes": {"pocket": [{"domain": r.domain_b, "sim": _f(r.kissim_sim, 3), "davis": _f(r.davis_coinhib, 3)}
                                  for r in look_pocket.itertuples()],
                       "sequence": nearest_seq[:6]},
        "range": {"n_chains": int(len(pdb)), "n_blind_esm": int(pdb.blind_esm.fillna(False).astype(bool).sum()) if len(pdb) else 0,
                  "n_novel_esm": int(pdb.novel_esm.fillna(False).astype(bool).sum()) if len(pdb) else 0,
                  "ligands": [{"ccd": c, "type": t, "n": n} for (c, t), n in lig_counts.most_common(12)],
                  "first_release": str(pdb.release_date.min()) if len(pdb) else None},
    }


# ------------------------------------------------------------------ driver
def nearest_by_sequence() -> dict[str, list]:
    p = DATA / "phylo" / "seq_tree.mldist"
    if not p.exists():
        return {}
    lines = p.read_text().split("\n")
    n = int(lines[0].strip())
    names, rows = [], []
    for line in lines[1:n + 1]:
        f = line.split()
        names.append(f[0])
        rows.append([float(x) for x in f[1:]])
    D = np.array(rows)
    out = {}
    for i, a in enumerate(names):
        order = np.argsort(D[i])
        out[a] = [{"domain": names[j], "dist": round(float(D[i, j]), 3)} for j in order[1:7]]
    return out


def run() -> dict:
    (APP / "data").mkdir(parents=True, exist_ok=True)
    JDIR.mkdir(parents=True, exist_ok=True)
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    s = pd.read_parquet(DATA / "annotation" / "structures.parquet")
    for col in ("tag", "construct_start", "is_rep", "ligand", "klifs_ligand_type", "pdb_id", "chain",
                "resolution", "release_date", "blind_af2", "blind_esm", "novel_af2", "novel_esm", "ptm", "sample_idx"):
        if col not in s:
            s[col] = None
    dom = pd.read_parquet(DATA / "annotation" / "domains.parquet").set_index("domain_id")
    motif = _read(DATA / "annotation" / "motif_dev.parquet")
    ribbon = _read(DATA / "annotation" / "ribbon.parquet")
    sites = _read(DATA / "pockets" / "sites.parquet")
    sim = _read(DATA / "selectivity" / "similarity.parquet")
    chains = _read(DATA / "pdb" / "chains.parquet")
    controls = _json(DATA / "pockets" / "controls.json", [])
    ifit = _read(DATA / "induced_fit" / "induced_fit.parquet")
    pci = _json(DATA / "annotation" / "pci.json", {})
    q5 = _json(DATA / "selectivity" / "q5.json", {})
    mr = _read(DATA / "annotation" / "motif_repair.parquet")
    nseq = nearest_by_sequence()

    # ---------------- structures: choose + write superposed crops
    if SDIR.exists():
        shutil.rmtree(SDIR)
    tasks = []
    keep_ids = {}
    for _, k in cat.iterrows():
        sd = s[s.domain_id == k.domain_id]
        ref = frame_for(k.to_dict(), sd)
        if ref is None:
            continue
        pdb_reps = sd[(sd.source == "pdb") & sd.is_rep.fillna(False).astype(bool)]
        pdb_reps = pdb_reps.assign(r=pd.to_numeric(pdb_reps.resolution, errors="coerce").fillna(9.9)).sort_values("r").head(12)
        es = sd[(sd.source == "esm2") & (sd.status == "done") & (sd.tag.fillna("apo") == "apo")]
        want = list(pdb_reps.struct_id) + list(sd[(sd.source == "afdb") & (sd.status == "done")].struct_id) + esm_keep(es)
        keep_ids[k.domain_id] = want
        rows = [r for r in sd[sd.struct_id.isin(want)].to_dict("records")]
        tasks.append((k.to_dict(), rows, ref))
    files: dict[str, dict] = {}
    with ProcessPoolExecutor(int(tuning().get("cpu_workers", 4))) as ex:
        for (kd, _, _), w in zip(tasks, ex.map(write_domain_structures, tasks, chunksize=4)):
            files[kd["domain_id"]] = w
    size_gb = sum(p.stat().st_size for p in SDIR.rglob("*.cif.gz")) / 1e9

    # ---------------- SQLite
    if DB.exists():
        DB.unlink()
    con = sqlite3.connect(DB)
    con.executescript(SCHEMA)
    kin_rows = [(k.domain_id, k.gene, k.uniprot, k["name"], k.group, k.family, int(k.domain_start), int(k.domain_end),
                 _b(k.is_dark), _b(k.is_pseudokinase), k.pseudo_evidence or None, int(k.n_pdb_chains),
                 int(dom.loc[k.domain_id, "n_pdb_states"]) if k.domain_id in dom.index else 0, k.priority_tier)
                for _, k in cat.iterrows()]
    con.executemany("INSERT INTO kinase VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", kin_rows)
    agree_map = {}
    for d_id, d in dom.iterrows():
        agree_map[f"afdb:{d_id}"] = d.get("afdb_agreement")
    st_rows = []
    for r in s[s.status == "done"].itertuples():
        path = files.get(r.domain_id, {}).get(r.struct_id)
        st_rows.append((r.struct_id, r.domain_id, r.source, r.pdb_id, r.chain,
                        None if pd.isna(r.sample_idx) else int(r.sample_idx), None if pd.isna(r.release_date) else str(r.release_date),
                        _f(r.resolution, 2), _b(r.blind_af2), _b(r.blind_esm), _b(r.novel_af2), _b(r.novel_esm),
                        r.kincore_spatial, r.kincore_dihedral, _b(r.ac_in), _f(r.d1, 2), _f(r.d2, 2), _f(r.mean_plddt, 2),
                        _f(r.ptm, 3), _f(r.pocket_rmsd_nearest, 3),
                        agree_map.get(r.struct_id) or (dom.loc[r.domain_id, "esm_agreement"]
                                                        if r.source == "esm2" and r.domain_id in dom.index
                                                        and r.struct_id == dom.loc[r.domain_id].get("esm_best") else None),
                        path, None if pd.isna(r.tag) else r.tag))
    con.executemany("INSERT INTO structure VALUES (" + ",".join("?" * 23) + ")", st_rows)
    ids = {r[0] for r in st_rows}
    if len(motif):
        m = motif[motif.struct_id.isin(ids)]
        con.executemany("INSERT INTO motif_dev VALUES (?,?,?,?)",
                        [(a, b, _f(c), _f(d_, 2)) for a, b, c, d_ in m[["struct_id", "motif", "rmsd", "mean_plddt"]].itertuples(index=False)])
    if len(ribbon):
        rb = ribbon[ribbon.struct_id.isin(ids)]
        con.executemany("INSERT INTO ribbon VALUES (?,?,?,?,?)",
                        [(a, int(b), None if pd.isna(c) else int(c), _f(d_, 2), _f(e, 1))
                         for a, b, c, d_, e in rb[["struct_id", "klifs_pos", "residue", "dev_A", "plddt"]].itertuples(index=False)])
    if len(chains):
        lig = chains[~chains.ligand.astype(str).isin(["0", "", "None"])]
        con.executemany("INSERT INTO ligand VALUES (?,?,?,?)",
                        [(f"pdb:{r.structure_ID}", str(r.ligand), r.klifs_ligand_type,
                          _b(any(bool(getattr(r, c, False)) for c in ("bp_I_A", "bp_I_B", "bp_II_in", "bp_II_out", "bp_II_B", "bp_III"))))
                         for r in lig.itertuples()])
    if len(sites):
        con.executemany("INSERT INTO site VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                        [(r.site_id, r.domain_id, r.source, _f(r.recovery), _f(r.p2rank_max), _f(r.vol_min, 1), _f(r.vol_max, 1),
                          r.subpocket, _b(r.known_in_pdb), _b(r.candidate), _f(r.rank_score, 4),
                          json.dumps({"residues": [int(x) for x in r.residues], "klifs": [int(x) for x in r.klifs_positions]}))
                         for r in sites.itertuples()])
    if len(sim):
        con.executemany("INSERT INTO similarity VALUES (?,?,?,?,?)",
                        [(a, b, c, _f(d_, 4), _f(e, 4)) for a, b, c, d_, e in
                         sim[["domain_a", "domain_b", "source", "kissim_sim", "davis_coinhib"]].itertuples(index=False)])
    for kind, name in (("sequence", "tree_sequence"), ("pocket_pdb", "tree_pocket_pdb"), ("pocket_esm", "tree_pocket_esm")):
        nwk = DATA / "phylo" / f"{name}.nwk"
        if nwk.exists():
            con.execute("INSERT INTO tree VALUES (?,?,?,?)", (name, kind, nwk.read_text(), f"static/data/{name}.json"))
    for c in controls:
        con.execute("INSERT INTO control VALUES (?,?,?,?,?,?)",
                    (c["gene"], c["domain_id"], c["site_name"],
                     json.dumps({"residues": c.get("residues"), "pdb": c.get("pdb_id"), "ligand": c.get("ligand"),
                                 "best_jaccard": c.get("best_jaccard"), "afdb_jaccard": c.get("afdb_jaccard"),
                                 "afdb_recovered": c.get("afdb_recovered"), "n_samples": c.get("n_samples"),
                                 "status": c.get("status")}),
                     _b(c.get("recovered")), _f(c.get("score"))))
    if len(ifit):
        done = ifit[ifit.status == "done"]
        con.executemany("INSERT INTO induced_fit VALUES (?,?,?,?,?,?)",
                        [(r.domain_id, r.ccd, r.pdb_id, _f(r.apo_frac_dfgout), _f(r.holo_frac_dfgout), _f(r.lig_rmsd, 2))
                         for r in done.itertuples()])

    # helper tables
    near_pocket = {}
    if len(sim):
        for src in ("pdb", "esm2"):
            for a, g in sim[sim.source == src].groupby("domain_a"):
                near_pocket.setdefault(a, g.sort_values("kissim_sim", ascending=False).iloc[0].domain_b)
    cand_count = sites[sites.candidate].groupby("domain_id").size().to_dict() if len(sites) else {}
    pdb_mode = s[(s.source == "pdb") & (s.status == "done") & (s.state != "Unassigned")].groupby("domain_id").state \
        .agg(lambda x: x.value_counts().index[0]).to_dict()
    tips = []
    for _, k in cat.iterrows():
        d = dom.loc[k.domain_id] if k.domain_id in dom.index else {}
        g = (lambda key: d.get(key) if hasattr(d, "get") else None)
        tips.append((k.domain_id, k.gene, k.uniprot, k.group, k.family, g("pdb_pip"), g("afdb_agreement"),
                     g("esm_agreement"), pdb_mode.get(k.domain_id), g("afdb_state"), g("esm_best_state"),
                     int(g("n_pdb_states") or 0), int(cand_count.get(k.domain_id, 0)), _b(k.is_dark), _b(k.is_pseudokinase),
                     (nseq.get(k.domain_id) or [{}])[0].get("domain"), near_pocket.get(k.domain_id)))
    con.executemany("INSERT INTO kinase_tips VALUES (" + ",".join("?" * 17) + ")", tips)
    for _, k in cat[cat.is_pseudokinase].iterrows():
        sub = mr[mr.domain_id == k.domain_id] if len(mr) else pd.DataFrame()
        call = lambda src: (sub[sub.source == src].call.value_counts().index[0] if len(sub) and (sub.source == src).any() else None)  # noqa: E731
        img = next(iter(files.get(k.domain_id, {}).values()), None)
        con.execute("INSERT INTO pseudo_view VALUES (?,?,?,?,?,?,?,?,?)",
                    (k.domain_id, k.gene, k.group, k.pseudo_evidence, call("afdb"), call("esm2"), call("pdb"),
                     json.dumps(sub.to_dict("records"), default=str) if len(sub) else "[]", img))
    for _, k in cat[cat.is_dark].iterrows():
        d = dom.loc[k.domain_id] if k.domain_id in dom.index else {}
        sd = s[(s.domain_id == k.domain_id) & (s.status == "done")]
        af = sd[sd.source == "afdb"]
        es = sd[(sd.source == "esm2") & (sd.tag.fillna("apo") == "apo")]
        ss = sites[(sites.domain_id == k.domain_id)].sort_values("rank_score", ascending=False) if len(sites) else pd.DataFrame()
        con.execute("INSERT INTO dark_view VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    (k.domain_id, k.gene, k.group, k.family, _f(d.get("esm_afdb_rmsd") if hasattr(d, "get") else None, 2),
                     _f(d.get("esm_entropy") if hasattr(d, "get") else None, 3),
                     d.get("esm_state_fracs") if hasattr(d, "get") else None,
                     _f(af.mean_plddt.iloc[0], 1) if len(af) else None, _f(es.mean_plddt.mean(), 1) if len(es) else None,
                     ss.iloc[0].site_id if len(ss) else None, int(ss.candidate.sum()) if len(ss) else 0))

    summary = {"pci": pci, "q5": q5, "build": {"date": pd.Timestamp.now().strftime("%Y-%m-%d"),
                                              "structures_gb": round(size_gb, 3)},
               "tanglegram": _json(DATA / "phylo" / "pocket_tree_info.json", {}),
               "sequence_tree": _json(DATA / "phylo" / "sequence_tree_info.json", {}),
               "findings": _json(DATA / "findings.json", {}), "references": _json(DATA / "references.json", []),
               "pipeline": _json(DATA / "pipeline.json", {}),
               "controls": controls,
               "induced_fit": ifit.to_dict("records") if len(ifit) else [],
               "motif_fidelity": _json(DATA / "annotation" / "motif_summary.json", {}),
               "construct_starts": {k.domain_id: int(max(1, k.domain_start - PAD)) for _, k in cat.iterrows()}}
    con.executemany("INSERT INTO summary VALUES (?,?)", [(k, json.dumps(v, default=str)) for k, v in summary.items()])
    con.commit()
    con.execute("VACUUM")
    con.close()

    # ---------------- specimen JSON
    family_of = dict(zip(cat.domain_id, cat.family))
    pdb_all = s[(s.source == "pdb") & (s.status == "done")]
    rng = np.random.default_rng(3)
    big = []
    for _, k in cat.iterrows():
        sd = s[s.domain_id == k.domain_id]
        fam_ids = [x for x, f in family_of.items() if f == k.family]
        fp = pdb_all[pdb_all.domain_id.isin(fam_ids)][["d1", "d2"]].dropna()
        if len(fp) > 500:
            fp = fp.iloc[rng.choice(len(fp), 500, replace=False)]
        fam_pts = [[round(float(a), 1), round(float(b), 1)] for a, b in fp.values]
        dd = {c: (None if isinstance(v, float) and np.isnan(v) else v)
              for c, v in (dom.loc[k.domain_id].to_dict() if k.domain_id in dom.index else {}).items()}
        if dd.get("pdb_states") is not None:
            dd["pdb_states"] = list(dd["pdb_states"])
        ks = sites[sites.domain_id == k.domain_id].sort_values("rank_score", ascending=False) if len(sites) else pd.DataFrame()
        site_list = [{"id": r.site_id, "recovery": _f(r.recovery), "p2rank": _f(r.p2rank_max), "subpocket": r.subpocket,
                      "known": bool(r.known_in_pdb), "candidate": bool(r.candidate),
                      "residues": [int(x) for x in r.residues]} for r in ks.itertuples()]
        simk = sim[sim.domain_a == k.domain_id].sort_values("kissim_sim", ascending=False) if len(sim) else pd.DataFrame(
            columns=["source", "domain_b", "kissim_sim", "davis_coinhib"])
        ch_dom = chains[chains.domain_id == k.domain_id] if len(chains) else pd.DataFrame(columns=["ligand", "klifs_ligand_type"])
        js = specimen_json(k.to_dict(), dd, sd, fam_pts, ribbon, site_list, simk, nseq.get(k.domain_id, []),
                           {sid: f"/structures/{p}" for sid, p in files.get(k.domain_id, {}).items()}, ch_dom)
        text = json.dumps(_clean(js), separators=(",", ":"), allow_nan=False, default=str)
        (JDIR / f"{k.domain_id}.json").write_text(text)
        if len(text) > 150_000:
            big.append((k.domain_id, len(text)))
    total_gb = (DB.stat().st_size + sum(p.stat().st_size for p in SDIR.rglob("*")) +
                sum(p.stat().st_size for p in (APP / "static" / "data").rglob("*"))) / 1e9
    log(f"bundle: SQLite {DB.stat().st_size / 1e6:.1f} MB, structures {size_gb:.2f} GB, total {total_gb:.2f} GB; "
        f"{len(big)} specimen JSON files over 150 KB")
    return {"total_gb": total_gb, "oversize_json": big}
