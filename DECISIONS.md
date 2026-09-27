# KINFOLK decisions

Each entry: decision, alternatives, reason, how to reverse.

## 2026-09-27 19:37  Domains absent from KLIFS
- Decision: Added 4 UniProt 'Protein kinase' domains with pockets mapped by hmmalign (PF00069) anchored on KLIFS-mapped domains
- Alternatives: Drop them; KLIFS-only catalogue
- Reason: Spec asks for >= 480 domains and the full human kinome
- How to reverse: Filter pocket_source == 'hmmalign' in kinases.parquet

## 2026-09-27 19:45  Catalytic residue calls for the pseudokinase flag
- Decision: a catalytic residue (beta3 Lys, HRD Asp, DFG Asp) counts as present if either the KLIFS pocket or the hmmalign (PF00069) alignment anchored to KLIFS numbering shows it. When both say it is absent, the evidence motif is taken from the HMM alignment where aligned, else from KLIFS. The KLIFS-only call is kept in `pseudo_evidence_klifs`.
- Alternatives: KLIFS pocket only; HMM only.
- Reason: KLIFS aligns the Tribbles xDFG one residue out of register (TRIB1 reads LES, the literature reads SLE); the HMM cannot align atypical PI3K/PI4K-like domains. The HMM misplaces the beta3 Lys in VRK1, VRK2 and BUB1, which KLIFS gets right. Requiring both sources to agree before flagging avoids either source's register errors.
- How to reverse: use `pseudo_evidence_klifs`, or drop `merge_calls` in catalogue.build.

## 2026-09-27 19:40  Domains absent from KLIFS
- Decision: Added 4 UniProt 'Protein kinase' domains with pockets mapped by hmmalign (PF00069) anchored on KLIFS-mapped domains
- Alternatives: Drop them; KLIFS-only catalogue
- Reason: Spec asks for >= 480 domains and the full human kinome
- How to reverse: Filter pocket_source == 'hmmalign' in kinases.parquet

## 2026-09-27 19:40  Wrong KLIFS UniProt accessions
- Decision: Reassign when the pocket matches another UniProt kinase at least 0.2 better (identity < 0.6 to the listed accession): GCK: P35557 (0.19) -> Q12851 (1.00)
- Alternatives: Keep the KLIFS accession; drop the entry
- Reason: The listed accession is a different protein
- How to reverse: Remove the reassignment block in catalogue.build

## 2026-09-27 19:40  Domains absent from KLIFS
- Decision: Added 4 UniProt 'Protein kinase' domains with pockets mapped by hmmalign (PF00069) anchored on KLIFS-mapped domains
- Alternatives: Drop them; KLIFS-only catalogue
- Reason: Spec asks for >= 480 domains and the full human kinome
- How to reverse: Filter pocket_source == 'hmmalign' in kinases.parquet

## 2026-09-27 19:40  Wrong KLIFS UniProt accessions
- Decision: Reassign when the pocket matches another UniProt kinase at least 0.2 better (identity < 0.6 to the listed accession): GCK: P35557 (0.19) -> Q12851 (1.00)
- Alternatives: Keep the KLIFS accession; drop the entry
- Reason: The listed accession is a different protein
- How to reverse: Remove the reassignment block in catalogue.build

## 2026-09-27 19:41  Domains absent from KLIFS
- Decision: Added 4 UniProt 'Protein kinase' domains with pockets mapped by hmmalign (PF00069) anchored on KLIFS-mapped domains
- Alternatives: Drop them; KLIFS-only catalogue
- Reason: Spec asks for >= 480 domains and the full human kinome
- How to reverse: Filter pocket_source == 'hmmalign' in kinases.parquet

## 2026-09-27 19:41  Wrong KLIFS UniProt accessions
- Decision: Reassign when the pocket matches another UniProt kinase at least 0.2 better (identity < 0.6 to the listed accession): GCK: P35557 (0.19) -> Q12851 (1.00)
- Alternatives: Keep the KLIFS accession; drop the entry
- Reason: The listed accession is a different protein
- How to reverse: Remove the reassignment block in catalogue.build

