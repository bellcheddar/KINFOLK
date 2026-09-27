"""hmmalign to the Pfam Pkinase HMM and anchor HMM columns to KLIFS positions.

Used (a) for domains that have no KLIFS pocket string and (b) for the
phylogeny MSA. Anchors are domains with a trusted KLIFS mapping: for every
KLIFS position we take the modal HMM match column across anchors.
"""
from __future__ import annotations

import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from . import DATA

HMM = DATA / "hmm" / "PF00069.hmm"


def hmmalign(seqs: dict[str, str], hmm: Path = HMM) -> dict[str, str]:
    """Return A2M rows: uppercase/'-' = match columns, lowercase/'.' = inserts."""
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "in.fa"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        out = subprocess.run(["hmmalign", "--outformat", "A2M", str(hmm), str(fa)],
                             capture_output=True, text=True, check=True).stdout
    rows: dict[str, list[str]] = {}
    name = None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            rows[name] = []
        elif name:
            rows[name].append(line.strip())
    return {k: "".join(v) for k, v in rows.items()}


def column_map(a2m: str, offset: int = 1) -> dict[int, int]:
    """HMM match column (1-based) -> residue number (index + offset)."""
    col, res, out = 0, 0, {}
    for c in a2m:
        if c.isupper():
            col += 1
            out[col] = res + offset
            res += 1
        elif c == "-":
            col += 1
        elif c.islower():
            res += 1
    return out


def anchor_columns(anchor_rows: dict[str, str], anchor_maps: dict[str, dict[int, int]],
                   offsets: dict[str, int]) -> dict[int, int]:
    """KLIFS position -> modal HMM column across anchor domains."""
    votes: dict[int, Counter] = {}
    for name, row in anchor_rows.items():
        res_to_col = {r: c for c, r in column_map(row, offsets[name]).items()}
        for pos, resnum in anchor_maps[name].items():
            if resnum in res_to_col:
                votes.setdefault(pos, Counter())[res_to_col[resnum]] += 1
    return {pos: c.most_common(1)[0][0] for pos, c in votes.items()}


def klifs_map_from_hmm(a2m: str, klifs_cols: dict[int, int], offset: int = 1) -> dict[int, int]:
    cmap = column_map(a2m, offset)
    return {pos: cmap[col] for pos, col in klifs_cols.items() if col in cmap}
