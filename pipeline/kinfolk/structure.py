"""Structure I/O and geometry helpers (gemmi based)."""
from __future__ import annotations


import gemmi
import numpy as np


def read_structure(path_or_text: str) -> gemmi.Structure:
    if str(path_or_text).lstrip().startswith("data_"):
        doc = gemmi.cif.read_string(str(path_or_text))
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        st = gemmi.read_structure(str(path_or_text))  # handles .gz for PDB and mmCIF
    st.setup_entities()
    st.remove_hydrogens()
    return st


def residues(st: gemmi.Structure, chain: str | None = None, model: int = 0):
    """Yield (chain_name, seqnum, gemmi.Residue) for polymer residues."""
    for ch in st[model]:
        if chain is not None and ch.name != chain:
            continue
        for res in ch:
            if res.het_flag == "H" and not gemmi.find_tabulated_residue(res.name):
                continue
            yield ch.name, res.seqid.num, res


def atom_table(st: gemmi.Structure, chain: str | None = None, atom: str = "CA") -> dict[int, np.ndarray]:
    out: dict[int, np.ndarray] = {}
    for _, num, res in residues(st, chain):
        a = res.find_atom(atom, "*")
        if a is not None and num not in out:
            out[num] = np.array(a.pos.tolist())
    return out


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (R, t) such that P @ R.T + t best fits Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


def rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    return float(np.sqrt(((P - Q) ** 2).sum(1).mean())) if len(P) else float("nan")


def superposed_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    R, t = kabsch(P, Q)
    return rmsd(P @ R.T + t, Q)


def dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return float(np.degrees(np.arctan2(y, x)))
