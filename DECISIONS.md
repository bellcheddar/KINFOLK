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

## 2026-09-27 19:47  AKT1 control construct
- Decision: Fold AKT1 PH + kinase (1-418, UniProt boundaries with 10-residue padding), 32 samples, in pass A
- Alternatives: Kinase domain only (cannot contain the interdomain site)
- Reason: Spec 12: the MK-2206 site lies between the PH and kinase domains
- How to reverse: Remove the phkd line from data/esm2/extra_jobs.jsonl

## 2026-09-27 19:50  IQ-TREE version
- Decision: IQ-TREE 3.1.3 (bioconda iqtree >= 2 resolved to 3.x), same options: -m LG+G4 -B 1000 -T AUTO --threads-max 16
- Alternatives: Pin iqtree 2.x
- Reason: IQ-TREE 3 is the maintained successor and accepts the spec options unchanged
- How to reverse: Set iqtree = "2.*" in pixi.toml

## 2026-09-27 20:03  Post-install fixes
- Decision: tools/post_install.sh after every environment change: unify libomp (torch vs conda) and clear the macOS hidden flag on .pth files
- Alternatives: KMP_DUPLICATE_LIB_OK=TRUE; non-editable install
- Reason: Duplicate OpenMP runtimes abort torch; Python 3.12.14 skips hidden .pth files so the editable kinfolk import vanished after the PyTDC solve
- How to reverse: Delete tools/post_install.sh and the post-install task

## 2026-09-27 20:03  Fingerprints: KiSSim-lite
- Decision: kissim cannot be installed (conda-forge kissim requires biopython <= 1.77, unavailable for Python 3.12; not on PyPI). Implemented KiSSim-lite per spec 7.2: 85 KLIFS positions x (size, H-bond donors, acceptors, charge, aromatic, aliphatic, side-chain exposure, distances to 4 KiSSim subpocket centres), labelled KiSSim-lite everywhere
- Alternatives: Separate Python 3.8 environment for kissim
- Reason: Spec tie-break rule; one environment keeps the pipeline simple
- How to reverse: Install kissim in a py38 pixi feature and swap selectivity.fingerprint

## 2026-09-27 20:12  Fold precision
- Decision: float32 on mps
- Alternatives: float32, bfloat16, float16
- Reason: Fastest dtype reproducing the float32 ABL1 pocket core to < 0.5 A: {"float32": {"seconds": 16.4, "core_rmsd": 0.0, "mean_plddt": 93.1, "ptm": 0.965}, "bfloat16": {"error": "process abort in MPSNDArrayMatrixMultiplication (destination and accumulator dtypes differ); observed 2026-09-27 20:10, esm 3.4.1.post1, torch 2.11.0"}, "float16": {"error": "not supported: ESMC precision options are bf16 and fp32 only, and trunk and ESMC must match"}}
- How to reverse: Edit fold_dtype in config/tuning.toml

## 2026-09-27 20:16  Fold settings
- Decision: num_loops 3, num_sampling_steps 100, float32 on MPS
- Alternatives: Library defaults (20 loops, 200 steps); 10 loops; 50 steps
- Reason: ABL1 against a crystal pocket: {"3/100": {"s": 16.6, "plddt": 93.1, "pocket85_rmsd_vs_crystal": 2.08}, "10/100": {"s": 35.4, "plddt": 93.6, "pocket85_rmsd_vs_crystal": 2.1}, "20/200": {"s": 67.1, "plddt": 93.5, "pocket85_rmsd_vs_crystal": 2.18}}. Defaults cost 4x with no pocket gain; 50 steps saves only ~10% (trunk dominates) and moves the pocket core 0.24 A
- How to reverse: Edit num_loops / num_sampling_steps in config/tuning.toml and clear data/manifests/fold.jsonl

## 2026-09-27 20:21  Sequence tree rooting
- Decision: Root on the edge maximising summed ePK group purity (members / MRCA leaves over AGC, CAMK, CK1, CMGC, STE, TK, TKL); outgroup side is 100% non-ePK
- Alternatives: Midpoint rooting (splits CMGC; purity 4.46/7); rooting on the atypical MRCA (atypicals are scattered long branches, not a clade)
- Reason: Spec 6.1 asks for the split that best separates ePK groups from atypicals; purity rose to 5.26/7
- How to reverse: Replace root_tree with root_at_midpoint in phylo.py

## 2026-09-27 20:24  GOBSMACKED deep link
- Decision: Render "Send pocket to GOBSMACKED" disabled with a tooltip
- Alternatives: Link anyway (would open GOBSMACKED without the pocket)
- Reason: gobsmacked.mdeller.com is live but none of its scripts (stages.js, viewer.js, clip.js, sequence_track.js) read query parameters, so ?uniprot=&pocket= is ignored (checked 2026-09-27)
- How to reverse: Change needed in GOBSMACKED: on load, read URLSearchParams uniprot (UniProt accession) and pocket (comma-separated UniProt residue numbers), fetch the AFDB model for the accession and pre-select those residues as the pocket. Then set GOBSMACKED_LINK_ENABLED = True in app/static/js/specimen.js

## 2026-09-27 20:34  Q6 fold order
- Decision: Pass Q6 (apo + holo, 8 samples each) runs after the minimum viable set and T2, before the extra T3 samples
- Alternatives: Strict tier order (T4 last)
- Reason: Q6 is a headline question; extra T3 samples only refine ensembles
- How to reverse: Move pass Q6 in foldqueue.build_jobs

## 2026-09-27 20:34  Q6 fold order
- Decision: Pass Q6 (apo + holo, 8 samples each) runs after the minimum viable set and T2, before the extra T3 samples
- Alternatives: Strict tier order (T4 last)
- Reason: Q6 is a headline question; extra T3 samples only refine ensembles
- How to reverse: Move pass Q6 in foldqueue.build_jobs

## 2026-09-27 20:49  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-27 20:52  Positive-control site resolution
- Decision: Resolve each control from KLIFS allosteric-ligand chains at run time: try ligands by frequency; reject crystallisation additives (PEG, glycols, buffers, ions) and ligands under 15 heavy atoms; use the ligand copy that contacts the kinase chain; apply a class rule (ABL1 and PDPK1 < 50% of site residues inside the 85 KLIFS positions, EGFR < 70%, AKT1 must touch the PH domain, MAP2K1 any). Chosen: ABL1 6NPE KWD, MAP2K1 4MNE 573, PDPK1 3ORZ 2A2, EGFR 4G5J 0WM, AKT1 4EJN 0R4
- Alternatives: Hand-typed residue lists (forbidden by spec 12); most frequent allosteric ligand without filters (picked PE5 and DTD additives)
- Reason: KLIFS allosteric annotations include additives; the class rules keep each site where its ligand class binds
- How to reverse: Edit ADDITIVES / SITE_RULE in pockets.py

## 2026-09-27 20:53  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-27 21:39  README polish
- Decision: README.md written to references/readme_template.md (header, badges, contact table, Why it matters intro, tables, author footer) and linted for its hard rules (no em dashes, rel=noopener on every target=_blank, no hype words)
- Alternatives: readme_forge.py
- Reason: readme_forge.py is referenced by the marcs-vibe-coding skill but is not present anywhere on this machine
- How to reverse: Run readme_forge.py over README.md when it is available

