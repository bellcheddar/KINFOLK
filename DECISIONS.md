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

## 2026-09-27 22:00  GOBSMACKED deep link (supersedes the disabled button)
- Decision: GOBSMACKED /prepare now reads ?uniprot=&pocket= (change in GOBSMACKED app/static/js/prepare.js: validates the accession and up to 200 residue numbers, clears the EGFR example and its forced 1M17 hint, fetches the AlphaFold DB model, prefills the pocket and computes the box). KINFOLK links to https://gobsmacked.mdeller.com/prepare with the full KLIFS pocket in UniProt numbering (pocket_uniprot in each specimen JSON), opening in a new tab
- Alternatives: Keep the button disabled; send only the four field-mark residues
- Reason: Requested by Marc; AFDB models use UniProt numbering, so KINFOLK residue numbers carry over unchanged. Tested in headless Chrome against the real prepare.html and prepare.js with stubbed APIs
- How to reverse: Revert the prepare.js change in GOBSMACKED and restore the disabled button in specimen.html/specimen.js. The live link needs GOBSMACKED redeployed

## 2026-09-27 22:00  README polish (supersedes the manual-only entry)
- Decision: Ran readme_forge.py (found at Vibe_Coding/marcs-vibe-coding/scripts/) with readme_meta.json; adopted README.polished.md as README.md; README.html, README.preview.html and elementor-widget-markdown.css generated for the site
- Alternatives: Template-only README
- Reason: The script exists locally, just not inside the installed skill folder
- How to reverse: Re-run readme_forge.py after FINDINGS.md is final if the README should quote headline numbers

## 2026-09-27 22:00  Claude Code permissions
- Decision: Wrote .claude/settings.json verbatim from spec Appendix A (scoped allow list; denies push, SSH, rsync, sudo, root/home rm -rf, playwright install)
- Alternatives: Leave permissions to the session defaults
- Reason: Requested by Marc; lets a tmux-launched session run unattended within scoped limits
- How to reverse: Delete .claude/settings.json

## 2026-09-27 22:07  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 00:12  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 01:34  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 03:36  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 05:38  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 07:39  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 09:40  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 11:40  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 12:03  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 12:05  Launching long jobs from zsh
- Decision: Start the fold queue and orchestrator with setopt NO_BG_NICE (or from tmux) so they run at nice 0
- Alternatives: Default zsh background launch (BG_NICE: nice +5)
- Reason: At nice +5 the fold queue crawled overnight while the Mac was idle (one batch 6.7 h instead of ~45 s), consistent with macOS deprioritising low-priority background work; at nice 0 it runs at full speed
- How to reverse: None needed; relaunch as documented in README Usage

## 2026-09-28 14:05  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 16:05  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 17:46  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 18:11  Fold worker recycling
- Decision: Supervisor runs fold workers in fresh processes, recycled every 40 batches; gc + torch.mps.empty_cache after each batch; an OOM at batch size 1 recycles the worker and retries once before failing
- Alternatives: One long-lived process (spec 2.2 fold_instances = 1)
- Reason: MPS kept per-shape graph allocations: other allocations grew to 156 GB over ~12 h and every batch then failed; a fresh process starts clean and costs ~15 s of model load per 40 batches
- How to reverse: Set RECYCLE_JOBS = 0 behaviour by running foldqueue.run() directly

## 2026-09-28 18:11  Import path
- Decision: PYTHONPATH=$PIXI_PROJECT_ROOT/pipeline set in pixi [activation.env]
- Alternatives: Rely on the editable .pth file
- Reason: macOS re-applied the hidden flag to the .pth (Python 3.12.14 skips hidden .pth files), silently breaking every spawned worker from 14:05; the environment lives in iCloud-synced Documents, which may be touching file flags
- How to reverse: Remove [activation.env] from pixi.toml

## 2026-09-28 18:11  Orchestrator guards
- Decision: Watchdog restarts the fold worker if no batch finishes for 30 min; the final phases run only if at most 2% of batches failed, otherwise a gate opens
- Alternatives: Liveness check only
- Reason: The overnight nice +5 crawl and the OOM cascade both went unnoticed by a liveness-only check
- How to reverse: Edit STALL_MIN / MAX_FAILED_FRAC in finish.py

## 2026-09-28 18:16  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 20:22  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-09-28 21:45  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-11-17 10:52  Tanglegram untangling
- Decision: Greedy child swaps on the pocket tree against the fixed sequence tree (3 passes, accept only swaps that reduce crossings)
- Alternatives: Rotate both trees; exact optimisation
- Reason: Spec 6.3 asks for a greedy heuristic; the sequence tree stays in its canonical order
- How to reverse: Change passes in phylo.untangle

## 2026-11-17 10:58  Pocket trees with incomplete distance matrices
- Decision: Drop the domain with the most NaN distances, repeatedly, until the matrix is complete, then run neighbour joining and record the dropped names in pocket_tree_info.json
- Alternatives: Impute missing distances with the matrix maximum; lower the 50-row threshold in fp_distance; skip the pocket trees entirely
- Reason: A NaN appears when two KiSSim-lite fingerprints share fewer than 50 KLIFS rows, and in practice it is caused by a few domains that map to almost no KLIFS positions (EIF2AK4_C for PDB; PLK5, STK40, TAF1_C, TAF1_N for ESMFold2). Dropping them leaves 303 and 518 leaves with real distances, where imputation would invent branch lengths
- How to reverse: Change kinfolk.phylo.drop_incomparable, or relax the row threshold in kinfolk.selectivity.fp_distance

## 2026-11-17 10:58  Tanglegram scroll region focusable
- Decision: Give the Plate III tanglegram container tabindex 0 so keyboard users can scroll it
- Alternatives: Shrink the figure to fit; move the scroll to the page
- Reason: The tanglegram is 303 rows tall inside a 760 px box, which axe flags as scrollable-region-focusable once the trees exist
- How to reverse: Remove the tabindex in app/templates/resemblance.html

## 2026-11-17 11:51  MASTL folded on a trimmed construct
- Decision: Fold residues 25 to 225 joined by GSGS to 735 to 845 instead of the full 821-residue construct
- Alternatives: Leave MASTL with the single sample it already had; fold on CPU (hours per sample); lower the sampling steps
- Reason: ESMFold2's SVD fails to converge on the full construct on MPS, twice, and MASTL's non-conserved middle region is disordered in AlphaFold DB (pLDDT below 70 from 221 to 740), so it contributes nothing to the pocket while costing most of the compute
- How to reverse: Delete data/esm2/models/MASTL/apo/s01 to s07 and refold the full construct on CPU

## 2026-11-17 11:53  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-11-17 11:54  Tanglegram untangling
- Decision: Greedy child swaps on the pocket tree against the fixed sequence tree (3 passes, accept only swaps that reduce crossings)
- Alternatives: Rotate both trees; exact optimisation
- Reason: Spec 6.3 asks for a greedy heuristic; the sequence tree stays in its canonical order
- How to reverse: Change passes in phylo.untangle

## 2026-11-17 11:55  Subpocket assignment
- Decision: Map site residues to KLIFS regions and label front / gate / BP-I / BP-II by majority; < 30% of residues inside the 85 KLIFS positions = outside KLIFS (allosteric); ATP cleft = >= 50% inside KLIFS
- Alternatives: KLIFS subpocket centroids (KiFrag)
- Reason: Region-based and reproducible from the KLIFS numbering already computed
- How to reverse: Edit SUBPOCKET in pockets.py

## 2026-11-17 11:55  Tanglegram untangling
- Decision: Greedy child swaps on the pocket tree against the fixed sequence tree (3 passes, accept only swaps that reduce crossings)
- Alternatives: Rotate both trees; exact optimisation
- Reason: Spec 6.3 asks for a greedy heuristic; the sequence tree stays in its canonical order
- How to reverse: Change passes in phylo.untangle

## 2026-11-17 11:57  Tanglegram untangling
- Decision: Greedy child swaps on the pocket tree against the fixed sequence tree (3 passes, accept only swaps that reduce crossings)
- Alternatives: Rotate both trees; exact optimisation
- Reason: Spec 6.3 asks for a greedy heuristic; the sequence tree stays in its canonical order
- How to reverse: Change passes in phylo.untangle

## 2026-11-17 11:58  PDB representative for pocket fingerprints
- Decision: Pick the chain whose KLIFS pocket identity is highest (in tenths), then the best resolution, then struct_id
- Alternatives: Best resolution alone, as before; require a minimum pocket identity and drop domains that fail it
- Reason: A KiSSim-lite fingerprint can only measure a pocket that is resolved. Resolution alone chose LRRK2 chains resolving 0.49 and 0.24 of the KLIFS pocket over a 2.9 A chain resolving all of it, and ties on resolution fell to whichever chain the parallel annotation wrote first, so the pocket tree changed between identical runs. Q5 PDB rho rises from 0.30 to 0.33 and the predicted-structure fraction falls from 1.04 to 0.96
- How to reverse: Restore the sort in kinfolk.selectivity.fingerprints to resolution alone

## 2026-11-17 11:58  Fold manifest records are deduplicated by id
- Decision: annotate keeps the last record per job id
- Alternatives: Rewrite the manifest in place when a batch is refolded
- Reason: The manifest is append-only, so a refolded batch had two done records and its samples were annotated twice, which broke the bundle on a unique constraint
- How to reverse: Revert the _last dict in kinfolk.annotate.run

