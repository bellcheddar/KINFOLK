"""KLIFS pocket numbering: region layout, REST client and sequence mapping."""
from __future__ import annotations

import numpy as np
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

API = "https://klifs.net/api"

# KLIFS 85-position region layout (van Linden et al. 2014), 1-based inclusive.
REGIONS: dict[str, tuple[int, int]] = {
    "I": (1, 3), "g.l": (4, 9), "II": (10, 13), "III": (14, 19), "aC": (20, 30),
    "b.l": (31, 37), "IV": (38, 41), "V": (42, 44), "GK": (45, 45), "hinge": (46, 48),
    "linker": (49, 52), "aD": (53, 59), "aE": (60, 64), "VI": (65, 67), "c.l": (68, 75),
    "VII": (76, 78), "VIII": (79, 79), "xDFG": (80, 83), "a.l": (84, 85),
}

# Stable pocket core used as the superposition frame (spec section 5.2).
CORE_REGIONS = ("I", "II", "III", "IV", "V", "GK", "hinge", "VI", "VII")
CORE_POSITIONS = [p for r in CORE_REGIONS for p in range(REGIONS[r][0], REGIONS[r][1] + 1)]

# Motifs reported for per-motif deviation (spec section 5.3).
MOTIFS = ("g.l", "aC", "b.l", "GK", "hinge", "c.l", "xDFG", "a.l")

# Landmark positions (KLIFS numbering).
POS_B3_LYS = 17
POS_AC_GLU = 24
POS_GATEKEEPER = 45
POS_HRD_ASP = 70     # the D of HRD sits at c.l position 70 in KLIFS numbering
POS_DFG_X = 80
POS_DFG_ASP = 81
POS_DFG_PHE = 82
POS_DFG_GLY = 83


def region_of(pos: int) -> str:
    for name, (a, b) in REGIONS.items():
        if a <= pos <= b:
            return name
    raise ValueError(pos)


@retry(stop=stop_after_attempt(4), wait=wait_exponential(min=2, max=30))
def get(endpoint: str, **params) -> list | dict:
    r = requests.get(f"{API}/{endpoint}", params=params, timeout=60)
    if r.status_code == 400 and "No" in r.text:  # KLIFS returns 400 for empty results
        return []
    r.raise_for_status()
    return r.json()


def map_pocket_to_sequence(pocket: str, sequence: str, offset: int = 1) -> dict[int, int]:
    """Map KLIFS positions (1..85) to residue numbers in `sequence`.

    `pocket` is the 85-character KLIFS pocket string ('_' or '-' for gaps).
    Returns {klifs_pos: residue_number} where residue_number = index + offset.
    Uses a semi-global alignment of the ungapped pocket residues against the
    sequence, so pocket runs that are contiguous in sequence stay contiguous
    while insertions between KLIFS positions are tolerated.
    """
    from biotite.sequence import ProteinSequence
    from biotite.sequence.align import SubstitutionMatrix, align_optimal

    positions = [i + 1 for i, c in enumerate(pocket) if c not in "_-"]
    residues = "".join(c for c in pocket if c not in "_-")
    clean = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in residues)
    target = "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in sequence)
    matrix = SubstitutionMatrix.std_protein_matrix()
    ali = align_optimal(ProteinSequence(clean), ProteinSequence(target), matrix,
                        gap_penalty=(-12, -1), terminal_penalty=False, max_number=1)[0]
    mapping: dict[int, int] = {}
    for i_p, i_s in ali.trace:
        if i_p >= 0 and i_s >= 0:
            mapping[positions[i_p]] = int(i_s) + offset
    return mapping


def pocket_identity(pocket: str, sequence: str, mapping: dict[int, int], offset: int = 1) -> float:
    """Fraction of the ungapped pocket residues reproduced by the mapping."""
    hits = sum(1 for p, r in mapping.items() if sequence[r - offset] == pocket[p - 1])
    return hits / max(1, sum(1 for c in pocket if c not in "_-"))


def core_coords(ca_by_resnum: dict[int, np.ndarray], mapping: dict[int, int],
                positions=CORE_POSITIONS) -> tuple[list[int], np.ndarray]:
    keep = [p for p in positions if p in mapping and mapping[p] in ca_by_resnum]
    return keep, np.array([ca_by_resnum[mapping[p]] for p in keep])
