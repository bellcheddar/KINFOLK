"""Kincore conformational labels (Modi and Dunbrack, PNAS 2019), keyed to KLIFS positions.

Rules and cluster centroids are ported verbatim from Kincore-standalone
(github.com/vivekmodi/Kincore-standalone, modules/spatial_label.py and
modules/dihedral_label.py; dihedral cutoff 0.45). Residues are located by KLIFS
position instead of Kincore's own HMM, so every source (PDB, AFDB, ESMFold2)
is labelled through the same KLIFS mapping.

The alphaC in/out call follows spec 5.4 (Lys NZ to Glu carboxylate < 4.0 A);
Kincore's own C-helix rule (Lys CB to Glu CB <= 10 A) is returned too.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import gemmi
import numpy as np

from . import klifs
from .structure import dihedral

CUTOFF = 0.45
POS_GLU4 = klifs.POS_AC_GLU + 4

# Side-chain atom standing in for Phe CZ when the DFG Phe is substituted (Kincore table).
PHE_ATOM = {"PHE": "CZ", "ARG": "CZ", "LEU": "CG", "PRO": "CG", "ASN": "CG", "MET": "CE", "SER": "OG",
            "HIS": "NE2", "VAL": "CB", "ALA": "CB", "TRP": "CZ3", "TYR": "OH", "GLY": "CA"}

DFGINTER = {"BABtrans": (-80.20, 128.22, -117.47, 23.76, -85.16, 133.21)}
DFGOUT = {"BBAminus": (-138.56, -176.12, -144.35, 103.66, -82.59, -9.03)}
DFGIN_MINUS = {"BLAminus": (-128.64, 178.67, 61.15, 81.21, -96.89, 20.53),
               "ABAminus": (-111.82, -7.64, -141.55, 148.01, -127.79, 23.32),
               "BLBminus": (-134.79, 175.48, 60.44, 65.35, -79.44, 145.34)}
DFGIN_PLUS = {"BLAplus": (-119.24, 167.71, 58.94, 34.08, -89.42, -8.54),
              "BLBplus": (-125.28, 172.53, 59.98, 32.92, -85.51, 145.28)}
DFGIN_TRANS = {"BLBtrans": (-106.16, 157.24, 69.37, 21.33, -61.73, 134.56)}

DIHEDRAL_LABELS = ("BLAminus", "BLAplus", "ABAminus", "BLBminus", "BLBplus", "BLBtrans", "BABtrans", "BBAminus")
SPATIAL_LABELS = ("DFGin", "DFGinter", "DFGout")


@dataclass
class KincoreResult:
    spatial: str = "Unassigned"
    dihedral: str = "Unassigned"
    d1: float | None = None          # alphaC Glu(+4) CA to DFG-Phe CZ
    d2: float | None = None          # beta3 Lys CA to DFG-Phe CZ
    lys_glu_salt: float | None = None
    lys_glu_cb: float | None = None
    ac_in: bool | None = None
    dihedral_dist: float | None = None
    angles: tuple | None = None

    def as_dict(self) -> dict:
        return asdict(self)


def spatial_label(d1: float | None, d2: float | None) -> str:
    if d1 is None or d2 is None:
        return "Unassigned"
    label = "Unassigned"
    if d1 <= 11 and d2 >= 11:
        label = "DFGin"
    if d1 > 11 and d2 <= 14:
        label = "DFGout"
    if d1 <= 11 and d2 <= 11:
        label = "DFGinter"
    return label


def _cos_dist(angles: tuple, centroid: tuple) -> float:
    return (2 / 6) * sum(1 - math.cos(math.radians(a - c)) for a, c in zip(angles, centroid))


def dihedral_label(spatial: str, angles: tuple | None, chi1: float | None) -> tuple[str, float | None]:
    if spatial == "Unassigned" or angles is None or chi1 is None or any(a is None for a in angles):
        return "Unassigned", None
    chi = chi1 % 360
    table = None
    if spatial == "DFGin":
        table = DFGIN_MINUS if 240 < chi <= 360 else DFGIN_PLUS if 0 < chi <= 120 else DFGIN_TRANS
    elif spatial == "DFGinter" and 120 < chi <= 240:
        table = DFGINTER
    elif spatial == "DFGout" and 240 < chi <= 360:
        table = DFGOUT
    if table is None:
        return "Unassigned", None
    best, best_d = "Unassigned", 999.0
    for name, cen in table.items():
        d = _cos_dist(angles, cen)
        if d <= best_d:
            best_d = d
            if d <= CUTOFF:
                best = name
    return best, round(best_d, 3)


def _atom(res: gemmi.Residue | None, name: str) -> np.ndarray | None:
    if res is None:
        return None
    a = res.find_atom(name, "*")
    return np.array(a.pos.tolist()) if a is not None else None


def _dist(a, b) -> float | None:
    return None if a is None or b is None else float(np.linalg.norm(a - b))


def _chi1_atom(res: gemmi.Residue) -> str | None:
    for n in ("CG", "SG", "CG1", "OG", "OG1"):
        if res.find_atom(n, "*") is not None:
            return n
    return None


def label(chain: gemmi.Chain, klifs_map: dict[int, int]) -> KincoreResult:
    """Label one chain. `klifs_map` maps KLIFS position -> author residue number in `chain`."""
    by_num: dict[int, gemmi.Residue] = {}
    ordered: list[gemmi.Residue] = []
    for res in chain:
        if res.seqid.icode.strip():
            continue
        by_num.setdefault(res.seqid.num, res)
        ordered.append(res)
    index = {id(r): i for i, r in enumerate(ordered)}

    def at(pos: int) -> gemmi.Residue | None:
        n = klifs_map.get(pos)
        return by_num.get(n) if n is not None else None

    lys, glu, glu4 = at(klifs.POS_B3_LYS), at(klifs.POS_AC_GLU), at(POS_GLU4)
    x, asp, phe = at(klifs.POS_DFG_X), at(klifs.POS_DFG_ASP), at(klifs.POS_DFG_PHE)
    out = KincoreResult()
    phe_atom = PHE_ATOM.get(phe.name) if phe is not None else None
    phe_xyz = _atom(phe, phe_atom) if phe_atom else None
    out.d1 = _dist(_atom(glu4, "CA"), phe_xyz)
    out.d2 = _dist(_atom(lys, "CA"), phe_xyz)
    out.spatial = spatial_label(out.d1, out.d2)

    nz = _atom(lys, "NZ")
    oe = [a for a in (_atom(glu, "OE1"), _atom(glu, "OE2")) if a is not None]
    if nz is not None and oe:
        out.lys_glu_salt = min(float(np.linalg.norm(nz - o)) for o in oe)
        out.ac_in = out.lys_glu_salt < 4.0
    out.lys_glu_cb = _dist(_atom(lys, "CB"), _atom(glu, "CB"))

    def neighbour(res, step):
        if res is None:
            return None
        i = index.get(id(res))
        if i is None or not (0 <= i + step < len(ordered)):
            return None
        nb = ordered[i + step]
        return nb if abs(nb.seqid.num - res.seqid.num) == 1 else None

    def phi_psi(res):
        if res is None:
            return None, None
        prev, nxt = neighbour(res, -1), neighbour(res, 1)
        n, ca, c = _atom(res, "N"), _atom(res, "CA"), _atom(res, "C")
        pc, nn = _atom(prev, "C"), _atom(nxt, "N")
        phi = dihedral(pc, n, ca, c) if all(v is not None for v in (pc, n, ca, c)) else None
        psi = dihedral(n, ca, c, nn) if all(v is not None for v in (n, ca, c, nn)) else None
        return phi, psi

    xp, xs = phi_psi(x)
    dp, ds = phi_psi(asp)
    fp, fs = phi_psi(phe)
    chi1 = None
    if phe is not None and (g := _chi1_atom(phe)):
        pts = [_atom(phe, n) for n in ("N", "CA", "CB", g)]
        if all(p is not None for p in pts):
            chi1 = dihedral(*pts)
    out.angles = (xp, xs, dp, ds, fp, fs, chi1)
    out.dihedral, out.dihedral_dist = dihedral_label(out.spatial, (xp, xs, dp, ds, fp, fs), chi1)
    return out
