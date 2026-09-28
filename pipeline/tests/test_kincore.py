"""Phase 3 acceptance: our Kincore port vs the Kincore database on 20 PDB chains (>= 18/20)."""
import json

import pandas as pd
import pytest

from kinfolk import DATA
from kinfolk.kincore import dihedral_label, label, spatial_label
from kinfolk.mapping import kincore_reference, map_chain
from kinfolk.structure import read_structure

RAW = DATA / "pdb" / "raw"


def test_spatial_rules():
    assert spatial_label(8.0, 13.0) == "DFGin"
    assert spatial_label(13.0, 9.0) == "DFGout"
    assert spatial_label(9.0, 9.0) == "DFGinter"
    assert spatial_label(13.0, 16.0) == "Unassigned"
    assert spatial_label(None, 5.0) == "Unassigned"


def test_dihedral_centroids_self_classify():
    assert dihedral_label("DFGin", (-128.64, 178.67, 61.15, 81.21, -96.89, 20.53), -70.9)[0] == "BLAminus"
    assert dihedral_label("DFGout", (-138.56, -176.12, -144.35, 103.66, -82.59, -9.03), -69.4)[0] == "BBAminus"
    assert dihedral_label("DFGinter", (-80.20, 128.22, -117.47, 23.76, -85.16, 133.21), 181.4)[0] == "BABtrans"
    assert dihedral_label("DFGin", (-119.24, 167.71, 58.94, 34.08, -89.42, -8.54), 55.6)[0] == "BLAplus"


def sample_chains(n=20, seed=7):
    ref = kincore_reference()
    ref = ref[ref.dihedral.isin(["BLAminus", "BLAplus", "ABAminus", "BLBminus", "BLBplus", "BLBtrans",
                                 "BABtrans", "BBAminus"])]
    st = pd.DataFrame(json.loads((DATA / "pdb" / "klifs_structures.json").read_text()))
    st = st[(st.species == "Human") & (st.alt.isin(["", "A"]))]
    st["pdb"] = st.pdb.str.upper()
    st = st[st.structure_ID.map(lambda s: (RAW / f"{s}.pdb.gz").exists())]
    m = st.merge(ref, on=["pdb", "chain"]).drop_duplicates(["pdb", "chain"])
    # Stratify: up to 3 per state, then fill.
    picks = pd.concat([g.sample(min(3, len(g)), random_state=seed) for _, g in m.groupby("dihedral")])
    if len(picks) < n:
        rest = m.drop(picks.index)
        picks = pd.concat([picks, rest.sample(min(n - len(picks), len(rest)), random_state=seed)])
    return picks.head(n)


@pytest.mark.skipif(not RAW.exists(), reason="PDB files not downloaded")
def test_kincore_matches_database():
    picks = sample_chains()
    assert len(picks) >= 20, f"only {len(picks)} labelled chains available"
    results = []
    for _, r in picks.iterrows():
        s = read_structure(str(RAW / f"{r.structure_ID}.pdb.gz"))
        ch = next((c for c in s[0] if c.name == r.chain), s[0][0])
        kmap, _ = map_chain(ch, r.pocket)
        k = label(ch, kmap)
        results.append((r.pdb + r.chain, r.spatial, r.dihedral, k.spatial, k.dihedral))
    hits = sum(1 for _, rs, rd, s, d in results if (rs, rd) == (s, d))
    report = "\n".join(f"{a}: ref {b}/{c} ours {d}/{e}" for a, b, c, d, e in results)
    (DATA / "kincore" / "unit_test_report.txt").write_text(f"{hits}/{len(results)} match\n{report}\n")
    assert hits >= 18, report
