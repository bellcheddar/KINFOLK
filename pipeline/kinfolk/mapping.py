"""Map structure chains to KLIFS positions (spec 5.1)."""
from __future__ import annotations

import gemmi
import pandas as pd

from . import DATA, klifs

THREE_TO_ONE = {k: gemmi.find_tabulated_residue(k).one_letter_code.upper()
                for k in ("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
                          "MSE SEP TPO PTR CSO CME KCX HIC MLY").split()}
THREE_TO_ONE.update({"MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C", "CME": "C",
                     "KCX": "K", "HIC": "H", "MLY": "K"})


def chain_sequence(chain: gemmi.Chain) -> tuple[str, list[int]]:
    seq, nums = [], []
    for res in chain:
        if res.seqid.icode.strip():
            continue
        one = THREE_TO_ONE.get(res.name)
        if one is None or res.find_atom("CA", "*") is None:
            continue
        seq.append(one)
        nums.append(res.seqid.num)
    return "".join(seq), nums


def map_chain(chain: gemmi.Chain, pocket: str) -> tuple[dict[int, int], float]:
    """KLIFS position -> author residue number, plus pocket identity."""
    seq, nums = chain_sequence(chain)
    if not seq:
        return {}, 0.0
    m = klifs.map_pocket_to_sequence(pocket, seq, offset=0)
    ident = klifs.pocket_identity(pocket, seq, m, offset=0)
    return {p: nums[i] for p, i in m.items()}, ident


def model_map(domain_seq: str, pocket: str, construct_start: int, domain_start: int,
              first_resnum: int = 1) -> dict[int, int]:
    """For predicted models numbered 1..N over a construct starting at `construct_start`."""
    m = klifs.map_pocket_to_sequence(pocket, domain_seq, offset=0)
    shift = domain_start - construct_start + first_resnum
    return {p: i + shift for p, i in m.items()}


def kincore_reference() -> pd.DataFrame:
    """Kincore database labels (PK_labels_PDB.fasta headers)."""
    p = DATA / "kincore" / "labels.parquet"
    if p.exists():
        return pd.read_parquet(p)
    rows = []
    for line in (DATA / "kincore" / "PK_labels_PDB.fasta").read_text().splitlines():
        if not line.startswith(">"):
            continue
        f = line[1:].split("\t")
        if f[0].startswith("AF-") or len(f) < 8:
            continue
        rows.append({"pdb": f[0][:4].upper(), "chain": f[0][4:], "gene": f[1], "spatial": f[4],
                     "dihedral": f[5], "ligands": f[6], "ligand_types": f[7]})
    df = pd.DataFrame(rows)
    df.to_parquet(p, index=False)
    return df


_CHAINS = None


def load_mapped(row, domain) -> tuple:
    """Open an annotated structure row and return (gemmi.Chain, {klifs_pos: resnum}).

    `row` is a structures.parquet row; `domain` a catalogue row (needs
    klifs_pocket_seq and full_sequence).
    """
    global _CHAINS
    from .structure import read_structure

    st = read_structure(str(DATA.parent / row["file_path"]))
    if row["source"] == "pdb":
        if _CHAINS is None:
            _CHAINS = pd.read_parquet(DATA / "pdb" / "chains.parquet", columns=["structure_ID", "pocket"]) \
                .set_index("structure_ID").pocket.to_dict()
        ch = next((c for c in st[0] if c.name == row["chain"]), st[0][0])
        kmap, _ = map_chain(ch, _CHAINS[int(str(row["struct_id"]).split(":")[1])])
        return ch, kmap
    uni = klifs.map_pocket_to_sequence(domain["klifs_pocket_seq"], domain["full_sequence"])
    ch = next(c for c in st[0] if len(c) > 20)
    if row["source"] == "esm2":
        cs = int(row["construct_start"])
        return ch, {p: u - cs + 1 for p, u in uni.items()}
    return ch, uni
