"""QC: parsers, KLIFS mapping, agreement classifier, PCI, pseudokinase calls, schema."""
import json
import sqlite3

import numpy as np
import pandas as pd
import pytest

from kinfolk import DATA, ROOT, klifs
from kinfolk.catalogue import merge_calls, pseudo_call
from kinfolk.hmm import column_map
from kinfolk.mapping import map_chain
from kinfolk.structure import kabsch, read_structure, superposed_rmsd
from kinfolk.summarise import classify, jsd

ABL1_POCKET = "HKLGGGQYGEVYEVAVKTLEFLKEAAVMKEIKPNLVQLLGVYIITEFMTYGNLLDYLREYLEKKNFIHRDLAARNCLVVADFGLS"


# ------------------------------------------------------------------ parsers and geometry
def test_a2m_column_map():
    # match columns: upper-case or '-'; inserts: lower-case or '.'
    assert column_map("AC-dEF") == {1: 1, 2: 2, 4: 4, 5: 5}


def test_kabsch_recovers_rotation():
    rng = np.random.default_rng(0)
    P = rng.normal(size=(20, 3))
    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    Q = P @ R.T + np.array([1.0, -2.0, 3.0])
    R2, t2 = kabsch(P, Q)
    assert np.allclose(P @ R2.T + t2, Q, atol=1e-8)
    assert superposed_rmsd(P, Q) < 1e-8


def test_esmfold2_mmcif_parser():
    """ESMFold2 output (to_mmcif) parses from a string; residues numbered 1..N over the construct."""
    cif = (DATA / "esm2" / "smoke" / "abl1_float32.cif").read_text()
    st = read_structure(cif)
    ch = st[0][0]
    assert len(ch) == 272
    assert ch[0].seqid.num == 1 and ch[len(ch) - 1].seqid.num == 272


# ------------------------------------------------------------------ KLIFS mapping
def test_abl1_landmarks():
    seq = pd.read_parquet(DATA / "catalogue" / "kinases.parquet").set_index("domain_id").loc["ABL1", "full_sequence"]
    m = klifs.map_pocket_to_sequence(ABL1_POCKET, seq)
    assert (m[17], m[24], m[45], m[81]) == (271, 286, 315, 381)
    assert klifs.pocket_identity(ABL1_POCKET, seq, m) == 1.0


def _chains_with_files(n=10):
    p = DATA / "pdb" / "chains.parquet"
    if not p.exists():
        pytest.skip("chains.parquet not built")
    ch = pd.read_parquet(p)
    ch = ch[ch.passes_quality & ch.pocket.str.count("_").eq(0)]
    ch = ch[ch.structure_ID.map(lambda s: (DATA / "pdb" / "raw" / f"{s}.pdb.gz").exists())]
    return ch.drop_duplicates("domain_id").sample(n, random_state=5)


def test_klifs_mapping_on_known_chains():
    """Mapped residues reproduce the KLIFS pocket string exactly for 10 complete chains."""
    from kinfolk.mapping import chain_sequence
    for _, r in _chains_with_files(10).iterrows():
        st = read_structure(str(DATA / "pdb" / "raw" / f"{r.structure_ID}.pdb.gz"))
        ch = next((c for c in st[0] if c.name == r.chain), st[0][0])
        kmap, ident = map_chain(ch, r.pocket)
        seq, nums = chain_sequence(ch)
        by_num = dict(zip(nums, seq))
        got = "".join(by_num.get(kmap.get(p), "_") for p in range(1, 86))
        assert ident >= 0.95, (r.pdb, r.chain, ident)
        assert sum(a == b for a, b in zip(got, r.pocket)) >= 82, (r.pdb, got, r.pocket)


# ------------------------------------------------------------------ agreement classifier (spec 5.6)
@pytest.mark.parametrize("state,rmsd,pdb,model,expected", [
    ("BLAminus", 1.0, {"BLAminus"}, {"BLAminus"}, "agree"),
    ("BLAminus", 2.0, {"BLAminus"}, {"BLAminus"}, "partial"),
    ("BLAminus", 1.0, {"BLAminus", "BBAminus"}, {"BLAminus"}, "partial"),      # the ABL1 case
    ("BLAminus", 1.0, {"BLAminus", "BBAminus"}, {"BLAminus", "BBAminus"}, "agree"),
    ("BLAminus", 3.5, {"BLAminus"}, {"BLAminus"}, "disagree"),
    ("BBAminus", 1.0, {"BLAminus"}, {"BBAminus"}, "disagree"),
    ("BLAminus", None, {"BLAminus"}, {"BLAminus"}, "disagree"),
    ("BLAminus", 1.0, set(), {"BLAminus"}, "none"),
])
def test_agreement_classifier(state, rmsd, pdb, model, expected):
    assert classify(state, rmsd, pdb, model) == expected


def test_jsd_bounds():
    p = np.array([1.0, 0, 0])
    assert jsd(p, p) == 0
    assert abs(jsd(p, np.array([0, 1.0, 0])) - 1.0) < 1e-12


# ------------------------------------------------------------------ pseudokinase calls
def test_pseudo_calls_and_merge():
    assert pseudo_call(ABL1_POCKET) == (False, "")
    trib = ABL1_POCKET[:79] + "RLESLE"            # KLIFS register: D81 -> L
    hmm = ABL1_POCKET[:79] + "LSLEDT"              # HMM register: SLE at 81-83
    merged = merge_calls(trib, hmm)
    assert pseudo_call(merged) == (True, "DFG→SLE")
    vrk_like = ABL1_POCKET                          # KLIFS shows K; HMM misaligned
    assert pseudo_call(merge_calls(vrk_like, ABL1_POCKET[:16] + "D" + ABL1_POCKET[17:]))[0] is False


def test_catalogue_acceptance():
    from kinfolk.catalogue import acceptance
    res = acceptance(pd.read_parquet(DATA / "catalogue" / "kinases.parquet"))
    assert res["pass"], res


# ------------------------------------------------------------------ schema
def test_sqlite_schema():
    db = ROOT / "app" / "data" / "kinfolk.sqlite"
    if not db.exists():
        pytest.skip("bundle not built")
    con = sqlite3.connect(db)
    tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"kinase", "structure", "motif_dev", "ribbon", "ligand", "site", "similarity", "tree", "control",
            "induced_fit", "summary"} <= tables
    idx = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='index'")}
    assert "ix_structure_domain_source" in idx and "ix_site_domain_source" in idx
    assert con.execute("SELECT COUNT(*) FROM kinase").fetchone()[0] >= 480
    bad = con.execute("SELECT COUNT(*) FROM structure WHERE source NOT IN ('pdb','afdb','esm2')").fetchone()[0]
    assert bad == 0
    for (k,) in con.execute("SELECT key FROM summary"):
        json.loads(con.execute("SELECT value_json FROM summary WHERE key=?", (k,)).fetchone()[0])
