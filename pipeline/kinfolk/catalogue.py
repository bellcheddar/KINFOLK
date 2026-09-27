"""Phase 1: human kinome catalogue, one row per catalytic domain (spec 4.1)."""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

from . import DATA, klifs
from . import hmm as hmmmod
from .config import Manifest, decision, log, tuning

RAW = DATA / "catalogue" / "raw"
OUT = DATA / "catalogue" / "kinases.parquet"
CONTROLS = {"ABL1", "MAP2K1", "PDPK1", "EGFR", "AKT1"}
KINASE_DOMAIN = re.compile(r"(protein kinase|kinase domain|PI3K/PI4K catalytic|alpha-type protein kinase|"
                           r"histidine kinase|ABC1|RIO kinase|FAT|Bromo)", re.I)
PSEUDO_TEXT = re.compile(r"pseudokinase|lacks? (the )?(protein )?kinase activity|catalytically (in)?active"
                         r"|devoid of (protein )?kinase activity|kinase-dead|lacks? (catalytic|key catalytic)", re.I)
JAK = {"JAK1", "JAK2", "JAK3", "TYK2"}
AA = set("ACDEFGHIKLMNPQRSTVWY")


@retry(stop=stop_after_attempt(4), wait=wait_exponential(min=2, max=30))
def _get_json(url: str, **params):
    r = requests.get(url, params=params, timeout=120)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def load_klifs() -> list[dict]:
    p = RAW / "klifs_kinase_information.json"
    if not p.exists():
        names = klifs.get("kinase_names", species="HUMAN")
        info = klifs.get("kinase_information", kinase_ID=",".join(str(x["kinase_ID"]) for x in names))
        p.write_text(json.dumps(info))
    return json.loads(p.read_text())


def load_uniprot(extra: set[str]) -> dict[str, dict]:
    p = RAW / "uniprot_kinases.json"
    entries = {e["primaryAccession"]: e for e in json.loads(p.read_text())}
    missing = sorted(extra - set(entries))
    px = RAW / "uniprot_extra.json"
    cached = {e["primaryAccession"]: e for e in json.loads(px.read_text())} if px.exists() else {}
    need = [a for a in missing if a not in cached]
    fields = "accession,gene_primary,protein_name,length,sequence,ft_domain,cc_similarity,cc_function,cc_domain,cc_caution,xref_hgnc"
    for i in range(0, len(need), 100):
        q = " OR ".join(f"accession:{a}" for a in need[i:i + 100])
        res = _get_json("https://rest.uniprot.org/uniprotkb/search", query=q, format="json",
                        fields=fields, size=500)
        for e in (res or {}).get("results", []):
            cached[e["primaryAccession"]] = e
    px.write_text(json.dumps(list(cached.values())))
    entries.update(cached)
    return entries


def _text(entry: dict) -> str:
    parts = []
    for c in entry.get("comments", []):
        for t in c.get("texts", []):
            parts.append(t.get("value", ""))
    return " ".join(parts)


def _gene(entry: dict) -> str:
    g = entry.get("genes") or [{}]
    return (g[0].get("geneName") or {}).get("value") or entry["primaryAccession"]


def _name(entry: dict) -> str:
    pd_ = entry.get("proteinDescription", {})
    return ((pd_.get("recommendedName") or {}).get("fullName") or {}).get("value", "")


def _domains(entry: dict) -> list[tuple[int, int, str]]:
    out = []
    for f in entry.get("features", []):
        if f["type"] == "Domain" and KINASE_DOMAIN.search(f.get("description", "")):
            out.append((f["location"]["start"]["value"], f["location"]["end"]["value"], f["description"]))
    return out


def _similarity_family(entry: dict) -> tuple[str, str]:
    for c in entry.get("comments", []):
        if c.get("commentType") == "SIMILARITY":
            t = " ".join(x["value"] for x in c.get("texts", []))
            fam = re.findall(r"([A-Za-z0-9/ -]+?) (?:sub)?family", t)
            return ("Other", fam[-1].strip() if fam else "Unclassified")
    return ("Other", "Unclassified")


def pseudo_call(pocket: str) -> tuple[bool, str]:
    """Objective pseudokinase flag from the KLIFS-aligned pocket (spec 4.1)."""
    def at(p):
        c = pocket[p - 1] if len(pocket) >= p else "-"
        return "-" if c in "_-" else c
    ev = []
    if at(klifs.POS_B3_LYS) != "K":
        ev.append(f"K→{at(klifs.POS_B3_LYS)}")
    if at(klifs.POS_HRD_ASP) != "D":
        ev.append("HRD→" + "".join(at(p) for p in (68, 69, 70)))
    if at(klifs.POS_DFG_ASP) != "D":
        ev.append("DFG→" + "".join(at(p) for p in (81, 82, 83)))
    return bool(ev), "; ".join(ev)


CHECKS = ((klifs.POS_B3_LYS, (klifs.POS_B3_LYS,), "K"),
          (klifs.POS_HRD_ASP, (68, 69, klifs.POS_HRD_ASP), "D"),
          (klifs.POS_DFG_ASP, (klifs.POS_DFG_ASP, 82, 83), "D"))


def merge_calls(klifs_pocket: str, hmm_pocket: str) -> str:
    """Per catalytic residue: canonical if either alignment shows it; otherwise
    report the HMM motif where aligned (literature register), else KLIFS."""
    out = list(klifs_pocket)
    for key, span, canon in CHECKS:
        k_ok, h_ok = klifs_pocket[key - 1] == canon, hmm_pocket[key - 1] == canon
        if k_ok:
            continue
        src = hmm_pocket if (h_ok or hmm_pocket[key - 1] not in "_-") else klifs_pocket
        for p in span:
            out[p - 1] = src[p - 1]
    return "".join(out)


def pdb_coverage(uniprots: list[str]) -> dict[str, list[dict]]:
    man = Manifest("pdbe_best_structures")
    cache_p = RAW / "pdbe_best_structures.json"
    cache = json.loads(cache_p.read_text()) if cache_p.exists() else {}
    todo = [u for u in uniprots if u not in cache]

    def fetch(u):
        try:
            d = _get_json(f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{u}")
            return u, (d or {}).get(u, [])
        except Exception as e:  # noqa: BLE001
            return u, {"error": str(e)[:200]}

    with ThreadPoolExecutor(tuning().get("io_workers", 8)) as ex:
        for u, rows in ex.map(fetch, todo):
            if isinstance(rows, dict):
                man.add(u, status=f"failed:{rows['error']}")
                continue
            cache[u] = rows
            man.add(u, n=len(rows))
    cache_p.write_text(json.dumps(cache))
    return cache


def domain_ids(rows: list[dict]) -> None:
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)
    for gene, rs in by_gene.items():
        rs.sort(key=lambda r: r["domain_start"])
        if len(rs) == 1:
            rs[0]["domain_id"] = gene
        elif gene in JAK and len(rs) == 2:
            rs[0]["domain_id"], rs[1]["domain_id"] = f"{gene}_JH2", f"{gene}_JH1"
        elif len(rs) == 2:
            rs[0]["domain_id"], rs[1]["domain_id"] = f"{gene}_N", f"{gene}_C"
        else:
            for i, r in enumerate(rs, 1):
                r["domain_id"] = f"{gene}_D{i}"


def build() -> pd.DataFrame:
    kl = [k for k in load_klifs() if k.get("pocket") and set(k["pocket"]) - set("_-")]
    uni = load_uniprot({k["uniprot"] for k in kl if k.get("uniprot")})
    log(f"catalogue: {len(kl)} KLIFS human entries with pockets, {len(uni)} UniProt entries")

    rows: list[dict] = []
    claimed: set[tuple[str, int, int]] = set()
    low_identity, skipped, reassigned = [], [], []
    for k in kl:
        e = uni.get(k["uniprot"])
        if e is None or "sequence" not in e:
            skipped.append(f"{k['name']} ({k.get('uniprot')})")
            continue
        seq = e["sequence"]["value"]
        m = klifs.map_pocket_to_sequence(k["pocket"], seq)
        ident = klifs.pocket_identity(k["pocket"], seq, m)
        if ident < 0.6:
            # KLIFS occasionally carries the wrong accession (e.g. GCK -> glucokinase).
            best = max(((klifs.pocket_identity(k["pocket"], x["sequence"]["value"],
                                               mm := klifs.map_pocket_to_sequence(k["pocket"], x["sequence"]["value"])), acc, mm)
                        for acc, x in uni.items() if "sequence" in x), key=lambda t: t[0])
            if best[0] >= ident + 0.2:
                reassigned.append(f"{k['name']}: {k['uniprot']} ({ident:.2f}) -> {best[1]} ({best[0]:.2f})")
                k = {**k, "uniprot": best[1]}
                e, seq, m, ident = uni[best[1]], uni[best[1]]["sequence"]["value"], best[2], best[0]
        if ident < 0.5:
            skipped.append(f"{k['name']} ({k['uniprot']}, pocket identity {ident:.2f}: wrong accession)")
            continue
        if ident < 0.8:
            low_identity.append((k["name"], round(ident, 2)))
        span = sorted(m.values())
        if len(span) < 20:
            skipped.append(f"{k['name']} (pocket maps only {len(span)} residues)")
            continue
        lo, hi = span[len(span) // 10], span[-len(span) // 10 - 1]
        doms = [d for d in _domains(e) if d[0] <= lo and hi <= d[1]]
        if doms:
            start, end = doms[0][0], doms[0][1]
        else:
            start, end = max(1, min(span) - 30), min(len(seq), max(span) + 45)
        claimed.add((k["uniprot"], start, end))
        rows.append({"gene": _gene(e), "uniprot": k["uniprot"], "name": _name(e), "klifs_name": k["name"],
                     "group": k["group"], "family": k["family"], "domain_start": start, "domain_end": end,
                     "sequence": seq[start - 1:end], "full_sequence": seq, "klifs_pocket_seq": k["pocket"],
                     "pocket_source": "klifs", "pocket_identity": round(ident, 3), "_entry": e})

    # UniProt kinase domains without a KLIFS entry: map via hmmalign anchored on KLIFS domains.
    extra = []
    for acc, e in uni.items():
        if "sequence" not in e:
            continue
        for s, t, desc in _domains(e):
            if not desc.lower().startswith("protein kinase"):
                continue
            if any(u == acc and not (t < a or s > b) for u, a, b in claimed):
                continue
            extra.append((acc, s, t, e))
    for a, s_, t, e in extra:
        grp, fam = _similarity_family(e)
        rows.append({"gene": _gene(e), "uniprot": a, "name": _name(e), "klifs_name": None,
                     "group": grp, "family": fam, "domain_start": s_, "domain_end": t,
                     "sequence": e["sequence"]["value"][s_ - 1:t], "full_sequence": e["sequence"]["value"],
                     "klifs_pocket_seq": None, "pocket_source": "hmmalign", "pocket_identity": None, "_entry": e})

    # Collapse duplicates (two KLIFS entries resolving to the same UniProt domain).
    best: dict[tuple, dict] = {}
    for r in rows:
        key = (r["uniprot"], r["domain_start"], r["domain_end"])
        if key not in best or (r["pocket_identity"] or 0) > (best[key]["pocket_identity"] or 0):
            best[key] = r
    if len(best) < len(rows):
        log(f"catalogue: collapsed {len(rows) - len(best)} duplicate domain rows")
    rows = list(best.values())

    # hmmalign every domain to PF00069, anchored to KLIFS numbering on well-mapped KLIFS domains.
    for i, r in enumerate(rows):
        r["_key"] = f"d{i}"
    anchors = [r for r in rows if r["pocket_source"] == "klifs" and r["pocket_identity"] >= 0.95
               and r["group"] != "Atypical"]
    a2m = hmmmod.hmmalign({r["_key"]: r["sequence"] for r in rows})
    amaps = {r["_key"]: {p: res - r["domain_start"] + 1 for p, res in
                         klifs.map_pocket_to_sequence(r["klifs_pocket_seq"], r["full_sequence"]).items()}
             for r in anchors}
    cols = hmmmod.anchor_columns({k: a2m[k] for k in amaps}, amaps, {k: 1 for k in amaps})
    for r in rows:
        m = hmmmod.klifs_map_from_hmm(a2m[r["_key"]], cols)
        r["hmm_pocket_seq"] = "".join(r["sequence"][m[p] - 1] if p in m else "_" for p in range(1, 86))
        if r["klifs_pocket_seq"] is None:
            r["klifs_pocket_seq"] = r["hmm_pocket_seq"]
    if extra:
        log(f"catalogue: {len(extra)} UniProt kinase domains without KLIFS entry mapped by hmmalign")
        decision("Domains absent from KLIFS", f"Added {len(extra)} UniProt 'Protein kinase' domains with "
                 "pockets mapped by hmmalign (PF00069) anchored on KLIFS-mapped domains",
                 "Drop them; KLIFS-only catalogue", "Spec asks for >= 480 domains and the full human kinome",
                 "Filter pocket_source == 'hmmalign' in kinases.parquet")
    if reassigned:
        log("catalogue: KLIFS accessions reassigned by pocket identity: " + "; ".join(reassigned))
        decision("Wrong KLIFS UniProt accessions", "Reassign when the pocket matches another UniProt kinase "
                 "at least 0.2 better (identity < 0.6 to the listed accession): " + "; ".join(reassigned),
                 "Keep the KLIFS accession; drop the entry", "The listed accession is a different protein",
                 "Remove the reassignment block in catalogue.build")
    if skipped:
        log(f"catalogue: skipped {len(skipped)} KLIFS entries with obsolete or missing UniProt records: "
            + ", ".join(skipped))
    if low_identity:
        log(f"catalogue: {len(low_identity)} KLIFS pockets map below 80% identity to UniProt canonical: "
            + ", ".join(f"{n} ({i})" for n, i in low_identity[:20]))

    domain_ids(rows)

    # Pseudokinase flags plus UniProt text cross-check.
    disagreements = []
    for r in rows:
        r["is_pseudokinase"], r["pseudo_evidence"] = pseudo_call(merge_calls(r["klifs_pocket_seq"], r["hmm_pocket_seq"]))
        r["pseudo_evidence_klifs"] = pseudo_call(r["klifs_pocket_seq"])[1]
        r["pseudo_uniprot_text"] = bool(PSEUDO_TEXT.search(_text(r["_entry"])))
        if r["is_pseudokinase"] != r["pseudo_uniprot_text"]:
            disagreements.append(f"{r['domain_id']} (objective={r['is_pseudokinase']}, "
                                 f"UniProt text={r['pseudo_uniprot_text']}, {r['pseudo_evidence'] or 'canonical'})")
    log(f"catalogue: pseudokinase flag disagreements with UniProt text: {len(disagreements)}")
    (DATA / "catalogue" / "pseudo_disagreements.txt").write_text("\n".join(disagreements) + "\n")

    # PDB coverage (SIFTS best structures via PDBe).
    cov = pdb_coverage(sorted({r["uniprot"] for r in rows}))
    for r in rows:
        dlen = r["domain_end"] - r["domain_start"] + 1
        chains = set()
        for c in cov.get(r["uniprot"], []):
            ov = min(c["unp_end"], r["domain_end"]) - max(c["unp_start"], r["domain_start"]) + 1
            if ov / dlen >= 0.7:
                chains.add((c["pdb_id"], c["chain_id"]))
        r["n_pdb_chains"] = len(chains)
        r["is_dark"] = len(chains) == 0
        r["n_pdb_states"] = None  # filled in Phase 3 (Kincore)

    for r in rows:
        if r["gene"] in CONTROLS:
            r["priority_tier"] = "T0"
        elif r["is_pseudokinase"] or r["is_dark"]:
            r["priority_tier"] = "T2"
        else:
            r["priority_tier"] = "T3"   # T1 (blind + novel) assigned in Phase 2
        r.pop("_entry"); r.pop("_key")

    df = pd.DataFrame(rows)
    cols = ["domain_id", "gene", "uniprot", "name", "group", "family", "domain_start", "domain_end",
            "sequence", "klifs_pocket_seq", "n_pdb_chains", "n_pdb_states", "is_dark", "is_pseudokinase",
            "pseudo_evidence", "priority_tier", "klifs_name", "pocket_source", "pocket_identity",
            "pseudo_uniprot_text", "pseudo_evidence_klifs", "hmm_pocket_seq", "full_sequence"]
    df = df[cols].sort_values("domain_id").reset_index(drop=True)
    dup = df["domain_id"].duplicated()
    if dup.any():
        df.loc[dup, "domain_id"] = df.loc[dup, "domain_id"] + "_" + df.loc[dup, "domain_start"].astype(str)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)
    log(f"catalogue: wrote {len(df)} domains; pseudokinases {int(df.is_pseudokinase.sum())}, "
        f"dark {int(df.is_dark.sum())}; groups {dict(Counter(df.group))}")
    return df


def acceptance(df: pd.DataFrame) -> dict:
    def ev(gene):
        s = df[df.gene == gene]
        return s.iloc[0]["pseudo_evidence"] if len(s) else "absent"
    res = {"n_domains": len(df), "n_ok": len(df) >= 480,
           "TRIB1": ev("TRIB1"), "ERBB3": ev("ERBB3"), "KSR2": ev("KSR2")}
    res["TRIB1_ok"] = "DFG→SLE" in res["TRIB1"]
    res["ERBB3_ok"] = "HRD→HRN" in res["ERBB3"]
    res["KSR2_ok"] = "K→R" in res["KSR2"]
    res["pass"] = all(res[k] for k in ("n_ok", "TRIB1_ok", "ERBB3_ok", "KSR2_ok"))
    return res
