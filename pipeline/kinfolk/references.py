"""About-page references (spec Appendix C): resolve every DOI via Crossref; drop and log failures."""
from __future__ import annotations

import json

import requests

from . import DATA
from .config import log

REFS = [
    {"key": "esmfold2", "what": "ESMFold2 / ESMC", "doi": "10.64898/2026.06.03.729735",
     "code": "https://github.com/Biohub/esm", "licence": "MIT"},
    {"key": "afdb", "what": "AlphaFold Protein Structure Database", "doi": "10.1093/nar/gkab1061",
     "code": "https://alphafold.ebi.ac.uk"},
    {"key": "klifs", "what": "KLIFS (2021 update)", "doi": "10.1093/nar/gkaa895", "code": "https://klifs.net"},
    {"key": "kincore", "what": "Kincore conformational classification", "doi": "10.1073/pnas.1814279116",
     "code": "https://github.com/vivekmodi/Kincore-standalone"},
    {"key": "kincore_db", "what": "Kincore database", "doi": "10.1093/nar/gkab920",
     "code": "http://dunbrack.fccc.edu/kincore"},
    {"key": "kissim", "what": "KiSSim (fingerprint design followed by KiSSim-lite)", "doi": "10.1021/acs.jcim.2c00050",
     "code": "https://github.com/volkamerlab/kissim"},
    {"key": "davis", "what": "Comprehensive kinase inhibitor profiling (Davis 2011)", "doi": "10.1038/nbt.1990"},
    {"key": "manning", "what": "The protein kinase complement of the human genome", "doi": "10.1126/science.1075762"},
    {"key": "p2rank", "what": "P2Rank", "doi": "10.1186/s13321-018-0285-8", "code": "https://github.com/rdk/p2rank"},
    {"key": "fpocket", "what": "fpocket", "doi": "10.1186/1471-2105-10-168", "code": "https://github.com/Discngine/fpocket"},
    {"key": "iqtree", "what": "IQ-TREE 2", "doi": "10.1093/molbev/msaa015", "code": "https://github.com/iqtree/iqtree3"},
    {"key": "mmseqs2", "what": "MMseqs2", "doi": "10.1038/nbt.3988", "code": "https://github.com/soedinglab/MMseqs2"},
    {"key": "molstar", "what": "Mol* Viewer", "doi": "10.1093/nar/gkab314", "code": "https://github.com/molstar/molstar"},
    {"key": "tdc", "what": "Therapeutics Data Commons (Davis loader)", "doi": "10.48550/arXiv.2102.09548",
     "code": "https://github.com/mims-harvard/TDC"},
]


def resolve() -> list[dict]:
    out, dropped = [], []
    for r in REFS:
        try:
            resp = requests.get(f"https://api.crossref.org/works/{r['doi']}", timeout=30,
                                headers={"User-Agent": "KINFOLK build (mailto:build@kinfolk.local)"})
            if resp.status_code != 200:
                raise ValueError(f"HTTP {resp.status_code}")
            m = resp.json()["message"]
            authors = m.get("author") or []
            first = authors[0].get("family", "") if authors else ""
            year = (m.get("issued") or {}).get("date-parts", [[None]])[0][0]
            out.append({**r, "title": (m.get("title") or [""])[0], "journal": (m.get("container-title") or [""])[0],
                        "first_author": first, "et_al": len(authors) > 2, "year": year})
        except Exception as e:  # noqa: BLE001
            dropped.append(f"{r['key']} ({r['doi']}: {e})")
    (DATA / "references.json").write_text(json.dumps(out, indent=1))
    log(f"references: {len(out)} DOIs resolved via Crossref; dropped {len(dropped)}: " + "; ".join(dropped))
    return out
