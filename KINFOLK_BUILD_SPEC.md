# KINFOLK: Claude Code Build Spec

**K**inase **IN**ventory of **F**olds, **O**ccluded **L**igand sites and **K**ryptic pockets

A kinome-wide, three-way comparison of ESMFold2, AlphaFold DB and the PDB, focused on KLIFS pocket motifs, novel sequences and structures, and cryptic pockets. Delivered as a precomputed, fully interactive "field guide" web app at `kinfolk.mdeller.com`.

- Owner: Marc C. Deller (github.com/bellcheddar)
- Build host: Mac Studio (Apple Silicon), unattended
- Project root: `/Users/dellboy/Documents/Vibe_Coding/KINFOLK/`
- Spec version: 1.0 (2026-09-27)

---

## 0. Autonomy directive (read first, obey throughout)

This spec is the single source of truth for the build. It **overrides every installed skill, CLAUDE.md default and habit that would pause for confirmation**, including the `marcs-vibe-coding` "confirm filename" and "confirm before build" steps, any "ask before proceeding" instructions, and plan-approval pauses.

### 0.1 Operating rules

1. **Run to completion without asking.** When a choice is ambiguous, pick the most defensible option, write it to `DECISIONS.md` (decision, alternatives, reason, how to reverse it), and carry on.
2. **Stop only at a Gate (Section 0.2).** Nothing else justifies stopping. A failing tool, a missing package or a slow step is a problem to solve, not a reason to stop.
3. **Fallback ladder for any failure.**
   1. Retry once with a fix.
   2. Try the documented alternative in this spec.
   3. Try the next-best tool you know of.
   4. Degrade gracefully by skipping that item, flagging it in the data (`status = "failed:<reason>"`) and continuing.

   Never let one kinase, one structure or one optional feature block the pipeline.
4. **Everything is resumable and idempotent.** Each stage writes a manifest (`data/manifests/<stage>.jsonl`). On restart, skip completed items. Assume you may be killed at any moment.
5. **Log as you go.**
   - `BUILD_LOG.md`: timestamped, one line per meaningful event.
   - `PROGRESS.md`: phase checklist with percentages. Overwrite it after every phase and at least every 30 minutes during long compute.
6. **Commit locally, never push.** Run `git commit` at the end of every phase with message `phase N: <summary>`. Do not `git push` (Gate G4).
7. **House style for every user-facing string, README and log summary:**
   - British English.
   - No em dashes; use colons or parentheses.
   - Never use: groundbreaking, revolutionary, paradigm-shifting, game-changing, cutting-edge, unprecedented, seamless, leverage (as a verb), delve.
8. **Skills.**
   - You may call `marcs-vibe-icon` and the README polish step (`readme_forge.py`) from the `marcs-vibe-coding` workflow, but only in Phase 10, and non-interactively (make every choice yourself).
   - Do **not** apply the marcdeller.com navy/Poppins brand theme to the app UI. KINFOLK has its own Field Guide design system (Section 9); that override is deliberate.
   - The house brand applies only to the README badge/icon.
9. **Stay inside the project root.** Never write outside `/Users/dellboy/Documents/Vibe_Coding/KINFOLK/`, apart from package caches (pixi, uv, pip, Hugging Face cache, Homebrew). Never use `sudo`.
10. **Keep the Mac awake.** Wrap long compute in `caffeinate -dimsu <command>`.

### 0.2 Gates (the only reasons to stop)

When a Gate opens:
1. Write `GATE_OPEN.md` with:
   - which gate;
   - what you need from Marc;
   - the exact command or answer that will unblock you;
   - what you have completed;
   - what you will do after the gate.
2. Send a macOS notification: `osascript -e 'display notification "<one line>" with title "KINFOLK gate <id>" sound name "Glass"'`.
3. **Keep working on everything the gate does not block**, and stop only when nothing unblocked remains.

| Gate | Opens when | Needs from Marc |
|------|------------|-----------------|
| **G1: Compute** | The ESMFold2 smoke test (Phase 0) fails on both MPS and CPU, or the projected fold time for the minimum viable set (Section 4.4) exceeds `MAX_COMPUTE_HOURS` (default 96) | Choose: accept the reduced scope, raise the budget, or supply a cloud GPU endpoint |
| **G2: Secrets** | Any credential is needed (Hugging Face token for a gated download, droplet SSH, GitHub) | The token/key, via environment variable, never pasted into files |
| **G3: Deploy** | Phase 10 is complete and the app passes QC locally | Approval to rsync to the droplet, plus the DNS/nginx/certbot steps |
| **G4: Publish** | Any `git push`, repo creation, or public release | Approval |
| **G5: Disk** | Free disk space on the project volume falls below 50 GB | Choose what to prune, or a new data volume |
| **G6: Science** | A positive control fails its acceptance criterion (Section 12) after the fallback ladder, **and** the failure would invalidate a headline figure | A decision on whether to proceed with the caveat flagged |

G6 exists to protect the headline claims, not to stop on noise. A single control failing, with the rest passing, is logged, not gated.

### 0.3 Launch recommendation (for Marc)

Before starting, copy Appendix A to `.claude/settings.json` and Appendix B to `CLAUDE.md` in the project root. Then launch inside tmux so the session survives a closed terminal:

```bash
cd /Users/dellboy/Documents/Vibe_Coding/KINFOLK
tmux new -s kinfolk
claude
# then: "Read KINFOLK_BUILD_SPEC.md and execute it end to end under Section 0."
```

The allowlist in Appendix A lets the build run unattended. The deny list blocks pushes, SSH, sudo and deletions outside the project, so autonomy stays scoped. Only run with all permission checks disabled if you are happy to accept that risk.

---

## 1. Scientific goals and headline questions

Every phase exists to answer one of these. Label every figure-level output with the question number.

- **Q1: Motif fidelity.** How accurately do ESMFold2 and AFDB reproduce each KLIFS motif (G-loop, β3 Lys, αC, gatekeeper, hinge, HRD, xDFG, activation loop) relative to the nearest experimentally observed state?
- **Q2: Conformational prior.** Do predictors collapse the kinome onto the active-like BLAminus state?
  - Quantify with the **Prior Collapse Index (PCI)**: the Jensen–Shannon divergence between the Kincore state distribution of predicted models and that of the PDB, per group and kinome-wide.
  - Test whether ESMFold2's diffusion ensemble recovers state diversity that a single AFDB model cannot.
- **Q3: Novel sequences and structures.** How do the models behave on:
  - the dark kinome (no PDB entry);
  - pseudokinases (the "motif repair" hallucination test);
  - post-cutoff blind structures (the hold-out sets)?
- **Q4: Cryptic pockets.** Do ESMFold2 ensembles recover known allosteric and cryptic sites from apo sequence alone (positive controls), and what new candidate sites emerge kinome-wide?
- **Q5: Selectivity.** Do pocket fingerprints computed from predicted structures predict inhibitor cross-reactivity as well as PDB-derived fingerprints do?
- **Q6: Induced fit (new capability).** ESMFold2 accepts ligands as CCD codes. When folded with a type II inhibitor, does the kinase adopt DFG-out? This is the causal test of the "apo prior" explanation behind Q2.

### Training cutoffs (define the blind sets)

| Model | Training PDB cutoff | Blind set = PDB entries released after |
|-------|--------------------|----------------------------------------|
| AlphaFold2 (AFDB) | 2018-04-30 | 2018-04-30 |
| ESMFold2 | September 2021 (model card: "Data Cutoff: Sept 2021") | 2021-09-30 |

A **sequence-novel** blind entry is one whose chain has < 40% identity (MMseqs2) to every chain released before the relevant cutoff. Both the "blind" and "blind + novel" tiers are reported.

---

## 2. Hardware optimisation (Mac Studio)

### 2.1 Phase 0 probe

Write `tools/hwprobe.py`. It records the following to `config/hardware.toml`:

| Item | Source |
|------|--------|
| Chip | `sysctl -n machdep.cpu.brand_string` |
| Performance cores | `sysctl -n hw.perflevel0.physicalcpu` |
| Efficiency cores | `sysctl -n hw.perflevel1.physicalcpu` |
| Unified memory | `sysctl -n hw.memsize` |
| GPU cores | `system_profiler SPDisplaysDataType` |
| Free disk | `df` |
| MPS availability | `torch.backends.mps.is_available()` |

Every downstream parallelism setting is derived from this file. **Never hard-code core or memory counts.**

### 2.2 Derived settings (`config/tuning.toml`, written by the probe)

| Setting | Rule |
|---------|------|
| `cpu_workers` | performance cores minus 2 (minimum 2) |
| `io_workers` | 16 (network-bound downloads, polite rate limits apply) |
| `fold_device` | `mps` if the smoke test passes, else `cpu` |
| `fold_dtype` | Try `bfloat16` on MPS first, then `float16`, then `float32`. Keep the fastest that reproduces the float32 reference ABL1 fold to < 0.5 Å Cα RMSD over the pocket core |
| `fold_instances` | 1 (7B parameters; one resident model, never reload per job) |
| `fold_batch_by_length` | Bucket sequences into 50-residue bins; process shortest first within each priority tier |
| `iqtree_threads` | `-T AUTO --threads-max <performance cores>` |
| `p2rank_threads` | `cpu_workers` |
| `memory_headroom_gb` | 24 (never let the fold process plus workers exceed total minus headroom) |

Set these environment variables for all PyTorch jobs:

```bash
export PYTORCH_ENABLE_MPS_FALLBACK=1      # unsupported ops fall back to CPU instead of crashing
export TOKENIZERS_PARALLELISM=false
export OMP_NUM_THREADS=<cpu_workers>
```

If memory pressure causes MPS allocation errors, set `PYTORCH_MPS_HIGH_WATERMARK_RATIO` to a tested value and log it in `DECISIONS.md`.

### 2.3 Concurrency architecture

The Studio's advantage is running the GPU fold queue and the CPU analysis pool **concurrently**:

```
          ┌──────────── fold queue (MPS, 1 resident ESMFold2) ─────────────┐
catalogue ┤                                                                ├─► annotation pool (CPU: Kincore, KLIFS mapping,
          └─ download pool (IO: PDB, AFDB, KLIFS) ─► prep pool (CPU) ──────┘    metrics, fpocket, P2Rank) ─► SQLite/Parquet
```

- Each fold writes its mmCIF and a `done` manifest line. The annotation pool polls the manifest and processes new models immediately, so analysis finishes minutes after the last fold rather than hours.
- Use `concurrent.futures.ProcessPoolExecutor` for CPU work and a single long-lived fold worker process.
- Communicate only through the filesystem and manifests: no shared memory, and no external queue service.
- Analytical store: Parquet plus DuckDB. App store: SQLite, read-only at serve time.

---

## 3. Environment

Two environments, matching the house pattern (pixi for compute as in BoltzMaker; uv for the Flask app as in GOBSMACKED).

### 3.1 Compute (`pixi.toml`, platform `osx-arm64`, channels `conda-forge`, `bioconda`)

| Package | Purpose | Fallback if unavailable on osx-arm64 |
|---------|---------|--------------------------------------|
| python 3.11 | runtime | none |
| pytorch (with MPS) | ESMFold2 | pip wheel |
| `esm` (pip; Biohub/esm, MIT) | ESMFold2 (`biohub/ESMFold2`), ESMC | Hugging Face `transformers` ≥ 5.16 ESMFold2 path |
| biotite, gemmi, numpy, scipy, pandas, pyarrow, duckdb | structure I/O, analysis | none |
| tmtools | TM-score | biotite superimposition |
| hmmer | hmmalign to Pfam Pkinase/PK_Tyr_Ser-Thr | Homebrew `hmmer` |
| mafft | MSA | Clustal Omega |
| iqtree (2.x) | ML phylogeny | FastTree |
| mmseqs2 | novelty and identity | Homebrew `mmseqs2`; else BLAST+ |
| fpocket | pocket detection | none (P2Rank covers it) |
| openjdk 17 + P2Rank release (github rdk/p2rank) | ML pocket prediction | fpocket only, flagged |
| opencadd | KLIFS remote client | direct KLIFS REST calls |
| kissim (volkamerlab) | KLIFS-pocket fingerprints | reimplement the published fingerprint (Section 7.2) |
| PyTDC | Davis kinase dataset | direct download of the Davis 2011 matrix |
| ete3 or scikit-bio | tree handling, NJ trees | Bio.Phylo |
| requests, tenacity, tqdm, rich | plumbing | none |

- **Kincore:** use Kincore-standalone (Modi & Dunbrack) if it installs cleanly. Otherwise implement the published rules in `kinfolk/kincore.py`:
  - spatial label (DFGin/DFGinter/DFGout) from the D1/D2 distances: DFG-Phe Cζ to the αC Glu(+4) Cα and to the β3 Lys Cα;
  - dihedral label (BLAminus, BLAplus, ABAminus, BLBminus, BLBplus, BLBtrans, BBAminus, BABtrans) from the X-D-F backbone and Phe χ1 dihedrals.

  Unit-test it against 20 KLIFS/Kincore-labelled PDB chains (Phase 3 acceptance).
- **PocketMiner** (TensorFlow) is **optional**. Attempt it only if `tensorflow-macos` installs in under 10 minutes; otherwise skip it and log the decision.

### 3.2 App (`app/pyproject.toml`, uv)

- Flask, gunicorn, Jinja2 and orjson. No build step and no Node toolchain.
- Front-end libraries are served as vendored files under `app/static/vendor/`, pinned and checksummed:
  - D3 v7;
  - Mol* viewer (molstar);
  - no other frameworks.

---

## 4. Data acquisition (Phases 1 and 2)

### 4.1 Kinome catalogue (Phase 1): `data/catalogue/kinases.parquet`

- **Sources:**
  - KLIFS (all human kinases, groups, families, 85-residue pocket sequences);
  - UniProt REST (canonical sequences; "Protein kinase" domain feature boundaries).
- **Unit of analysis = catalytic domain**, not gene. Multi-domain kinases get one row per domain (for example `JAK2_JH1`, `JAK2_JH2`).
- **Columns:** `domain_id, gene, uniprot, group, family, domain_start, domain_end, sequence, klifs_pocket_seq, n_pdb_chains, n_pdb_states, is_dark, is_pseudokinase, pseudo_evidence, priority_tier`.
- **`is_pseudokinase`** is computed objectively from the aligned sequence: flag the domain if any of these canonical catalytic residues is absent:
  - β3 Lys (VAIK);
  - catalytic Asp (HRD);
  - DFG Asp.

  Record which residue is missing in `pseudo_evidence` (for example "HRD→HRN", "DFG→SLE", "K→R"). Cross-check against UniProt/literature annotation where available and log disagreements.
- **`is_dark`:** no PDB chain covering at least 70% of the domain.

### 4.2 PDB (Phase 2a)

- **Source:** KLIFS processed structures (pocket-aligned, with DFG/αC annotations and ligand info), supplemented by RCSB for release dates and resolution.
- **Quality filter:** X-ray ≤ 3.0 Å or cryo-EM ≤ 3.5 Å, and pocket completeness ≥ 80% of the 85 KLIFS residues.
- **Stored per chain:**
  - Kincore labels;
  - αC in/out;
  - ligand CCD codes and KLIFS ligand type;
  - release date;
  - the `blind_af2`, `blind_esm`, `novel_af2` and `novel_esm` flags.
- **Deduplicate** to per-kinase state representatives: cluster pocket Cα coordinates at 1.0 Å RMSD and keep the highest-resolution member. Keep all chains for distributions.

### 4.3 AlphaFold DB (Phase 2b)

- **Source:** `https://alphafold.ebi.ac.uk/api/prediction/{uniprot}`. Download the model mmCIF and the PAE JSON.
- **Handling:** crop to the domain (keep the full-length file), and store per-residue pLDDT.
- Missing entries are flagged as `status="absent"`, not failed.

### 4.4 ESMFold2 (Phase 2c): the long pole

- **Model:** `EsmFold2Model.from_pretrained("biohub/ESMFold2", device=<fold_device>)`, 7B parameters, MIT licence.
- **API:** follow the API in the installed `esm` package; the reference pattern is `ESMFold2InputBuilder().fold(model, StructurePredictionInput(...), num_loops=…, num_sampling_steps=…, num_diffusion_samples=…)`.
- **Input:** domain sequence padded by 10 residues on each side (clipped to the protein ends).
- **Outputs per sample:** mmCIF (`to_mmcif()`), per-residue pLDDT, PAE and pTM.
- **Use Fold-CP** only for any domain construct above the length where MPS memory fails (expected to be rare for domain-only folds).
- **Ensembles come free:** ESMFold2 is a diffusion model, so `num_diffusion_samples` gives the conformational ensemble natively. No masking or dropout hacks are needed.

**Priority tiers** (fold in this order, so a partial run still yields the headline figures):

| Tier | Contents | Samples per domain |
|------|----------|-------------------|
| T0 | Positive controls (Section 12) + ABL1 | 32 |
| T1 | Blind + novel hold-out domains | 16 |
| T2 | Pseudokinases + dark kinome | 16 |
| T3 | All remaining domains | 8 |
| T4 | Q6 ligand-conditioned folds (Section 8) | 8 per complex |

**Minimum viable set for Gate G1:** T0 + T1 + one sample for every other domain.

**Budgeting:**
1. The Phase 0 smoke test times ABL1 (about 270 residues) at the default settings.
2. Extrapolate the total using length-squared scaling per bucket.
3. If the total exceeds `MAX_COMPUTE_HOURS`, reduce `num_sampling_steps` first (test 100 → 50 and check pocket-core RMSD against the 100-step fold stays < 0.5 Å), then reduce T3 samples down to 1.
4. Only if the minimum viable set still does not fit, open G1.

---

## 5. Annotation and metrics (Phase 3)

For every structure (PDB chain, AFDB model, each ESMFold2 sample):

1. **KLIFS mapping.** Align the domain to KLIFS pocket numbering (hmmalign to the kinase HMM, then map to the 85 KLIFS positions using KLIFS's own alignment for anchoring). Store `klifs_pos → residue` maps.
2. **Superposition frame.** Superpose on the **stable pocket core**: KLIFS regions I, II, III, IV, V, GK, hinge, VI and VII, Cα only. All motif deviations are measured in this frame, so a moving αC or activation loop cannot corrupt the fit.
3. **Per-motif deviation** relative to the **nearest PDB state representative** (the minimum over representatives, not a single reference). Motifs are g.l, αC, b.l, GK, hinge, c.l, xDFG and a.l, plus per-residue Cα deviation for all 85 positions (the pocket ribbon).
4. **Kincore labels** (spatial and dihedral), plus αC in/out by the Lys–Glu salt-bridge distance (< 4.0 Å in).
5. **Confidence:** mean pLDDT per motif; PAE between the N-lobe and C-lobe.
6. **Agreement class** per model (drives the tree pips):
   - **agree:** the model's Kincore label is in the set of labels observed in the PDB for that kinase **and** the pocket-core-frame 85-residue Cα RMSD to the nearest representative is < 1.5 Å.
   - **partial:** the label matches but the RMSD is 1.5–3.0 Å, **or** the model matches one state while the PDB records others the model never samples (the ABL1 case).
   - **disagree:** otherwise.
   - **none:** the kinase has no PDB structure (hollow pip).
   - **PDB pip** colour encodes state coverage instead: ≥ 2 Kincore states observed = full, 1 = partial, 0 = hollow.
7. **Ensemble summaries** (ESMFold2): the fraction of samples per Kincore state, state entropy, and whether the ensemble reaches any PDB state the single AFDB model misses ("ensemble rescue").
8. **Prior Collapse Index:** the Jensen–Shannon divergence between the pooled state distributions (PDB vs AFDB; PDB vs ESMFold2 single best sample; PDB vs ESMFold2 ensemble), kinome-wide and per group, with bootstrap 95% CIs over kinases (1,000 resamples).
9. **Pseudokinase "motif repair" test:** for each pseudokinase, does the model place the degraded residue in canonical catalytic geometry? For example, does a residue at the HRD Asp position sit where a catalytic Asp would, or does SLE sit in DFG-in geometry? Score it by the Cα/Cβ position error versus the canonical position in the closest active homologue. Output: restored / partial / faithful.

---

## 6. Phylogeny (Phase 4)

1. **Sequence tree.**
   - Build the MSA of all kinase domains with `hmmalign`; trim columns with > 50% gaps (log the count).
   - Run IQ-TREE 2: `-m LG+G4 -B 1000 -T AUTO`. If the run would take more than 6 hours, fall back to `-m LG+G4 -fast`, then to FastTree `-lg -gamma`.
   - Midpoint root, then reroot on the split that best separates the eukaryotic protein kinase groups from the atypicals; log the choice.
   - Output Newick with branch lengths and ultrafast bootstrap support.
2. **Pocket tree.** Compute pairwise KiSSim distances between kinase pocket fingerprints (PDB-representative based, and separately ESMFold2-best-sample based), then build neighbour-joining trees.
3. **Tanglegram.** Pair the sequence tree with the pocket tree. Minimise crossings by rotating internal nodes with a greedy heuristic, and record the crossing count before and after.
4. **Export for the front end:** `app/static/data/tree_sequence.json`, `tree_pocket_pdb.json`, `tree_pocket_esm.json`. Each is a nested JSON of `{name, length, support, children}` plus a leaf index into the kinase table. Precompute layout-independent data only; D3 does the layout (Section 9.3).

---

## 7. Cryptic pockets and selectivity (Phases 5 and 6)

### 7.1 Pocket detection (Q4)

- **Detection:** run fpocket and P2Rank on every ESMFold2 sample, the AFDB model and the PDB representatives. Use PocketMiner too if it was installed.
- **Clustering:** cluster pockets across an ensemble by residue-set Jaccard similarity (≥ 0.5, same site).
- **Per site, record:**
  - recovery frequency (the fraction of samples in which it appears);
  - maximum P2Rank probability;
  - volume range;
  - lining residues mapped to KLIFS subpockets (front, gate area, back pocket BP-I/BP-II) or "outside KLIFS" (allosteric);
  - whether the site is present in any PDB structure of that kinase ("known") or not ("candidate cryptic").
- **Candidate cryptic sites** need: recovery ≥ 0.25, absent from every PDB apo structure, P2Rank probability ≥ 0.5, and not in the ATP cleft. Rank candidates by recovery × probability × (1 + Open Targets tractability flag, where available).

### 7.2 Selectivity validation (Q5)

- **Fingerprints:** KiSSim fingerprints per kinase from three sources: PDB representative, AFDB, ESMFold2 best sample.
- **Ground truth:** the Davis 2011 profiling matrix (via PyTDC), converted to per-inhibitor binary hit calls at Kd ≤ 1 µM, giving kinase–kinase co-inhibition similarity (Jaccard over shared inhibitor hits).
- **Metric:** Spearman correlation between fingerprint similarity and co-inhibition similarity for each source, with bootstrap CIs. **Headline: the predicted-structure correlation as a fraction of the PDB correlation.**
- **Tie-break rule:** if kissim will not install, implement a simplified 85-position physicochemical fingerprint (residue type class, size, side-chain exposure, subpocket distances) and label it "KiSSim-lite" everywhere.

---

## 8. Induced-fit test (Phase 7, Q6)

1. From the PDB set, select up to 30 complexes that are **blind for ESMFold2** (released after 2021-09-30), DFG-out, with a type II ligand per the KLIFS ligand classification, and a ligand that has a CCD code. Prefer diversity across kinases.
2. Fold each twice: apo and with the ligand as CCD input (`LigandInput`; follow the installed API), with 8 samples each.
3. Readouts:
   - the fraction of samples in DFG-out, apo vs holo;
   - ligand heavy-atom RMSD to the crystal (after pocket-core superposition);
   - the Kincore label shift.
4. Include ABL1 with imatinib as the anchor case only if a post-cutoff ABL1/imatinib entry exists; otherwise report it separately as "memorisation possible".

---

## 9. Web app (Phases 8 and 9)

### 9.1 Architecture

- **Compute and serving are separate.** The app is a **precomputed atlas**: the Flask routes serve SQLite queries and static JSON, and no folding or heavy compute runs at request time. It must run on the 3.8 GB droplet.
- **Data bundle, `app/data/`:**
  - `kinfolk.sqlite` (tables in Section 10);
  - per-domain JSON under `static/data/kinase/<domain_id>.json`;
  - domain-cropped structures as gzipped mmCIF: the PDB representatives, the AFDB model, and at most 8 ESMFold2 samples per domain (the medoid of each Kincore state present, plus the best-pLDDT sample).
- **Size budget:** total bundle ≤ 2.5 GB. If it is over, reduce the stored samples to 4 and log the change.
- **URL state everywhere:** `?k=ABL1_JH&tree=phylo&src=seq&tips=agree&plate=2`. Every view is shareable and restorable.

### 9.2 Design system: Field Guide (from the approved B2 mockup)

**Colour tokens** (CSS custom properties in `static/css/tokens.css`):

| Token | Value | Use |
|-------|-------|-----|
| `--ground` | `#E6EADF` | page background |
| `--card` | `#F8F7F1` | plates and cards |
| `--ink` | `#1F3326` | text, sidebar, tree branches |
| `--ink-muted` | `#52625A` | captions |
| `--rule` | `#C9D0C1` | borders |
| `--sidebar-text` | `#F1EEDF` | sidebar text |
| `--agree` | `#3F6B52` | agreement pips |
| `--partial` | `#C98F2B` | partial, AFDB marker |
| `--disagree` | `#A9482A` | disagreement, selected-path highlight |
| `--esm` | `#E07A55` | ESMFold2 marker |
| `--band` | `#E9EDE2` | clade bands |
| `--band-selected` | `#F2E3DA` | selected clade band |
| `--ribbon-ok` | `#D9DFD0` | ribbon cell < 1 Å |

**Typography** (Google Fonts):
- Instrument Serif (display and headings, italic for taxonomy);
- Public Sans (body and UI);
- JetBrains Mono (identifiers, residue numbers, axis labels).

**Layout:**
- A fixed 288 px dark sidebar containing the wordmark, backronym, kinase search, the Plates navigation (Pl. I to VII) and the key to markers.
- The main column holds plates as cards with 4 px radius.

**Rules:**
- Colour is never the only carrier of meaning. Markers also differ in shape: PDB ring, AFDB square, ESMFold2 triangle.
- Contrast is ≥ 4.5:1 for text.
- Touch targets are ≥ 44 px.
- Every control is a real `<button>`, `<select>` or `<a>`; charts get `aria-label`s.
- Keyboard navigation works across tree tips (arrow keys move between leaves; Enter selects).

### 9.3 Plates (tabs)

| Plate | Route | Content |
|-------|-------|---------|
| **I. Family Tree** (home) | `/` | The B2 phylogram plus the selected specimen docked below |
| **II. Specimen Card** | `/k/<domain_id>` | Full per-kinase page |
| **III. Family Resemblance** | `/resemblance` | Tanglegram and off-target finder |
| **IV. Hidden Rooms** | `/rooms` | Cryptic pocket atlas |
| **V. Black Sheep** | `/black-sheep` | Pseudokinases |
| **VI. Distant Relatives** | `/dark` | Dark kinome |
| **VII. The Almanac** | `/almanac` | Scrollytelling findings |
| About | `/about` | Pipeline and references |

**Plate I: Family Tree** (build this first; it is the centrepiece)
- **Top-down phylogram in D3** (`d3.hierarchy` + `d3.cluster` for the cladogram; phylogram y-coordinates from cumulative branch length):
  - tips aligned with dotted leaders;
  - three pips per tip (PDB coverage, AFDB agreement, ESMFold2 agreement);
  - soft group bands;
  - a scale bar in substitutions per site.
- **Controls:**
  - Phylogram / Cladogram / Radial (radial uses `d3.cluster` with polar projection; animate between layouts with 600 ms transitions);
  - "Tree from": Sequence / KLIFS pocket (KiSSim) / Tanglegram;
  - "Tips show": Agreement / Kincore state / Cryptic pocket count / Dark kinome.
- **Scale:** about 500 tips will not fit legibly at a fixed width, so:
  - the default view collapses each group to a triangle sized by member count;
  - clicking a group expands it;
  - include a zoom/pan (`d3.zoom`) and a "fit group" button;
  - labels render only above a zoom threshold; otherwise use a hover tooltip.
- **Selection:** highlights the root-to-tip path in `--disagree` at 2.6 px and opens a tooltip card (nearest kin by sequence, nearest kin by pocket, pip summary, "Open specimen" link).
- **Bootstrap support:** optional small node labels for support < 70 (toggle).

**Plate II: Specimen Card**
- Header with gene, full name, taxonomy (Group › Family › UniProt › domain range) and actions ("Compare with a relative", "Send pocket to GOBSMACKED").
  - The GOBSMACKED action deep-links to `https://gobsmacked.mdeller.com/?uniprot=<id>&pocket=<klifs_residue_list>`. If that parameter contract does not yet exist, render the button disabled with a tooltip and log the change needed in GOBSMACKED in `DECISIONS.md`.
- **Conformational compass:** a D1/D2 scatter of every PDB chain in the kinase family (grey), with this kinase's PDB states (rings), the AFDB model (square) and all ESMFold2 samples (triangles, opacity by pLDDT). Cluster labels for DFGin, DFGinter and DFGout.
- **Field marks:** β3 Lys, αC Glu, gatekeeper and DFG Asp, with residue numbers (UniProt numbering).
- **Observer's note:** an auto-generated, template-based summary sentence (no LLM), for example "Both predictors settle this kinase in BLAminus; the PDB also records DFG-out, seen only with type II ligands."
- **KLIFS ribbon:** 85 cells × rows (AFDB, ESMFold2 best, ESMFold2 ensemble mean), tinted by deviation (< 1 / 1–2.5 / > 2.5 Å), with region labels.
- **Mol\* panel:** a three-way overlay (PDB representative, AFDB, ESMFold2) in the pocket-core frame. Hovering a ribbon cell highlights that residue in Mol\*, and clicking a residue in Mol\* highlights the ribbon cell (two-way linking via the Mol\* plugin API).
- **Cards:** Hidden rooms, Look-alikes (pocket kin vs sequence kin) and Range (PDB chains, ligands, blind entries).

**Plate III: Family Resemblance**
- A tanglegram of the sequence tree vs the pocket tree, with crossing lines coloured by group and discordant pairs emphasised.
- An **Off-target finder:** pick a kinase, and it lists kinases ranked by pocket similarity, per source (PDB / AFDB / ESMFold2), with the Davis co-inhibition score beside each.
- The **Q5 summary chart:** correlation per source with CIs.

**Plate IV: Hidden Rooms**
- A kinome-wide table of candidate cryptic sites, sortable and filterable.
- A **site viewer:** Mol\* with pocket blobs rendered as translucent surfaces, and a recovery-frequency threshold slider.
- A **KLIFS subpocket sunburst.**
- A **positive-control scorecard** (Section 12), always visible at the top.

**Plate V: Black Sheep**
- A pseudokinase gallery: cards that flip between the degraded motif (sequence) and the model's geometry (restored / partial / faithful), each with a Mol\* snapshot.

**Plate VI: Distant Relatives**
- The dark kinome: the domain list with ESMFold2–AFDB agreement as an orthogonal confidence score, ensemble state entropy, and top candidate pockets.

**Plate VII: The Almanac**
- Scrollytelling, with one section per question Q1 to Q6 and each figure pulled from live data.
- The **PCI figure:** a "slider morphs the kinome" animation of the Kincore state distribution, PDB → AFDB → ESMFold2 single → ESMFold2 ensemble.

**About**
- A pipeline schematic (SVG); a software and reference table with GitHub links and DOIs; data cutoffs; the licence; and how to cite.

### 9.4 Performance targets

- First contentful paint < 1.5 s on the droplet.
- Tree interactions at 60 fps with 500 tips.
- Specimen JSON ≤ 150 KB, gzipped via nginx.

---

## 10. Data model (SQLite, `app/data/kinfolk.sqlite`)

```sql
kinase(domain_id PK, gene, uniprot, name, grp, family, dom_start, dom_end,
       is_dark, is_pseudo, pseudo_evidence, n_pdb_chains, n_pdb_states, tier)
structure(struct_id PK, domain_id FK, source CHECK(source IN ('pdb','afdb','esm2')),
          pdb_id, chain, sample_idx, release_date, resolution,
          blind_af2, blind_esm, novel_af2, novel_esm,
          kincore_spatial, kincore_dihedral, ac_in, d1, d2,
          mean_plddt, ptm, pocket_rmsd_nearest, agreement, file_path)
motif_dev(struct_id FK, motif, rmsd, mean_plddt)
ribbon(struct_id FK, klifs_pos, residue, dev_A, plddt)
ligand(struct_id FK, ccd, klifs_type, in_backpocket)
site(site_id PK, domain_id FK, source, recovery, p2rank_max, vol_min, vol_max,
     subpocket, known_in_pdb, candidate, rank_score, residues_json)
similarity(domain_a, domain_b, source, kissim_sim, davis_coinhib)
tree(tree_id PK, kind, newick, json_path)
control(control_id PK, domain_id, site_name, expected_residues_json, recovered, score)
induced_fit(domain_id, ccd, pdb_id, apo_frac_dfgout, holo_frac_dfgout, lig_rmsd)
summary(key PK, value_json)   -- PCI values, correlations, headline numbers
```

Create indices on every foreign key and on `(domain_id, source)`.

---

## 11. Phase plan (execute in order; the fold queue runs in the background from Phase 2 onward)

| Phase | Name | Main outputs | Acceptance criteria |
|-------|------|--------------|---------------------|
| 0 | Bootstrap and probe | repo skeleton, envs, `hardware.toml`, `tuning.toml`, ESMFold2 ABL1 smoke test, fold-time projection | ABL1 folds; the dtype decision is logged; the projection is written; G1 checked |
| 1 | Catalogue | `kinases.parquet` | ≥ 480 human domains; the pseudokinase flags reproduce the known cases TRIB1 (DFG→SLE), ERBB3 (HRD→HRN) and KSR2 (β3 K→R) |
| 2 | Acquisition | PDB chains, AFDB models; **fold queue started** (tiers T0→T4) | ≥ 95% of domains have an AFDB model or an `absent` flag; blind flags populated |
| 3 | Annotation | Kincore, KLIFS maps, motif deviations, agreement, ensembles, PCI | Kincore unit test ≥ 18/20 match; ABL1 PDB shows both DFGin and DFGout states |
| 4 | Phylogeny | sequence tree, pocket trees, tanglegram JSON | Groups monophyletic except known exceptions (log them); runtime logged |
| 5 | Pockets | sites table, candidate list, control scorecard | Section 12 criteria evaluated; G6 checked |
| 6 | Selectivity | similarity table, Q5 correlations | The PDB-source Spearman ρ is positive with a CI excluding 0 (sanity check of the method) |
| 7 | Induced fit | `induced_fit` table | ≥ 10 complexes processed, or all available if fewer |
| 8 | Data bundle | SQLite, JSON, trimmed mmCIF | Bundle ≤ 2.5 GB; the schema validates; every kinase JSON loads |
| 9 | Web app | Flask app, all plates | The Section 13 QC passes locally |
| 10 | Polish and handover | README (readme_forge), icon (marcs-vibe-icon), About page, `FINDINGS.md`, deploy script, **Gate G3** | Everything committed; `GATE_OPEN.md` written for G3 |

**Build order inside Phase 9:**
1. Tokens and layout shell.
2. Plate I tree.
3. Plate II Specimen Card.
4. Plates IV, III, V, VI, VII, About.

Screenshot each plate after it is built (Section 13) before moving on.

**While waiting on folds:** if the fold queue is still running when a phase needs ESMFold2 data, build that phase against AFDB and PDB data first, mark its ESMFold2 fields "pending", and backfill automatically when the manifests show completion. Never idle while the queue is running.

---

## 12. Positive controls (Q4 acceptance)

Resolve each control's site residues from the cited ligand-bound PDB structures at run time (the residues contacting the ligand within 4 Å). **Do not hand-type residue lists.**

| Control | Site | Ligand class used to define the site |
|---------|------|---------------------------------------|
| ABL1 | Myristoyl pocket (C-lobe) | asciminib / GNF-2 class allosteric inhibitors |
| MEK1 (MAP2K1) | Type III pocket adjacent to ATP | allosteric MEK inhibitors (for example the trametinib/selumetinib class) |
| PDK1 (PDPK1) | PIF pocket (N-lobe) | PIF-pocket ligands |
| EGFR | EAI045 allosteric site | EAI045 |
| AKT1 | PH–kinase interdomain pocket | allosteric AKT inhibitors (MK-2206 class) |

**Criterion:** a control counts as recovered if an ESMFold2 site cluster from the **apo, sequence-only** ensemble overlaps the defined site with Jaccard ≥ 0.4 in at least 1 sample. Report the recovery frequency too.

**Note:** AKT1's interdomain site needs the PH domain, which the domain-only fold does not include. For this control, fold the PH + kinase construct (UniProt boundaries) as a special case, and log it in `DECISIONS.md`.

Report the same metric for AFDB (single model) as the baseline.

---

## 13. QC (Phase 9 and 10)

1. **Tests.** `pytest`:
   - parsers, the Kincore rules, KLIFS mapping on 10 known chains, and the agreement classifier on hand-built fixtures;
   - the SQLite schema, and route smoke tests (every route returns 200 for 20 random kinases).
2. **Visual.** Playwright (Chromium; the project's Playwright finds the system browser; do not run `playwright install`) screenshots each plate at 1440×900 and 390×844 into `qc/screens/`. Check each for:
   - no horizontal scroll;
   - no overlapping tree labels at default zoom;
   - no console errors;
   - the Mol\* viewer initialising.
3. **Accessibility.** axe-core via Playwright on each plate, with zero serious or critical violations.
4. **Science sanity checks,** written to `qc/science.md`:
   - The ABL1 story reproduces: PDB DFGin and DFGout; AFDB and ESMFold2 best sample are DFGin.
   - The pseudokinase flags match Phase 1.
   - The PCI has a CI.
   - Q5 PDB ρ > 0.
5. **Performance.** Lighthouse (or Playwright timing) FCP < 1.5 s for the local gunicorn serving the production bundle.

---

## 14. Deliverables and handover

In the project root at the end:

- `app/`: the deployable Flask app with its data bundle.
- `pipeline/`: all compute code (`kinfolk/` package, CLI entry points `kinfolk run --phase N`, `kinfolk status`).
- `FINDINGS.md`: headline results for Q1 to Q6, with numbers, CIs and the figure paths. Plain summary style, British English, no hype words.
- `README.md`: polished with readme_forge, with the icon and badge.
- `DECISIONS.md`, `BUILD_LOG.md`, `PROGRESS.md`, `qc/`.
- `deploy/`:
  - `deploy.sh` (rsync bundle → droplet, systemd unit `kinfolk.service` for gunicorn on 127.0.0.1:8010, nginx server block for `kinfolk.mdeller.com`, certbot command);
  - **Not executed.** Deployment waits on Gate G3.
- `GATE_OPEN.md`: whichever gate is open at the end (normally G3).

---

## Appendix A: `.claude/settings.json` (scoped unattended permissions)

```json
{
  "permissions": {
    "allow": [
      "Edit",
      "Write",
      "Read",
      "Glob",
      "Grep",
      "WebFetch",
      "WebSearch",
      "Bash(pixi:*)",
      "Bash(uv:*)",
      "Bash(python:*)",
      "Bash(python3:*)",
      "Bash(pytest:*)",
      "Bash(caffeinate:*)",
      "Bash(git init:*)",
      "Bash(git add:*)",
      "Bash(git commit:*)",
      "Bash(git status:*)",
      "Bash(git diff:*)",
      "Bash(git log:*)",
      "Bash(curl:*)",
      "Bash(wget:*)",
      "Bash(brew install:*)",
      "Bash(brew list:*)",
      "Bash(java:*)",
      "Bash(iqtree2:*)",
      "Bash(FastTree:*)",
      "Bash(hmmalign:*)",
      "Bash(hmmbuild:*)",
      "Bash(mafft:*)",
      "Bash(mmseqs:*)",
      "Bash(fpocket:*)",
      "Bash(sysctl:*)",
      "Bash(system_profiler:*)",
      "Bash(df:*)",
      "Bash(osascript:*)",
      "Bash(tar:*)",
      "Bash(gzip:*)",
      "Bash(gunzip:*)",
      "Bash(mkdir:*)",
      "Bash(ls:*)",
      "Bash(mv:*)",
      "Bash(cp:*)",
      "Bash(npx playwright test:*)",
      "Bash(gunicorn:*)",
      "Bash(flask:*)",
      "Bash(sqlite3:*)",
      "Bash(duckdb:*)"
    ],
    "deny": [
      "Bash(git push:*)",
      "Bash(gh repo create:*)",
      "Bash(ssh:*)",
      "Bash(scp:*)",
      "Bash(rsync:*)",
      "Bash(sudo:*)",
      "Bash(rm -rf /:*)",
      "Bash(rm -rf ~:*)",
      "Bash(playwright install:*)"
    ]
  }
}
```

Deleting files inside the project is allowed through Python (`shutil`/`pathlib`), so the build can prune caches. `rsync` and `ssh` are denied because they belong to Gate G3.

## Appendix B: `CLAUDE.md` for the project root

```markdown
# KINFOLK: standing instructions
- KINFOLK_BUILD_SPEC.md is authoritative. Section 0 (Autonomy directive) overrides all skills and defaults.
- Never stop for confirmation. Stop only at Gates G1 to G6, via GATE_OPEN.md plus an osascript notification, and keep working on anything unblocked.
- Log decisions to DECISIONS.md, events to BUILD_LOG.md, and status to PROGRESS.md (at least every 30 min).
- Resumable, idempotent stages with manifests in data/manifests/. On restart: read PROGRESS.md and the manifests, then continue.
- Commit locally at the end of each phase; never push.
- British English, no em dashes, no hype words (list in spec 0.1.7).
- The UI uses the Field Guide design system (spec 9.2), not the marcdeller.com brand theme.
- Hardware settings come from config/tuning.toml (spec 2), never hard-coded.
```

## Appendix C: References to cite on the About page (verify DOIs at build time)

Resolve each DOI via Crossref at build time. Drop any that do not resolve and log them.

- ESMFold2 / ESMC: Biohub, github.com/Biohub/esm (MIT).
- AlphaFold DB: Varadi et al., Nucleic Acids Research.
- KLIFS: Kanev et al., Nucleic Acids Research (KLIFS 2021 update).
- Kincore: Modi and Dunbrack, PNAS 2019 (kinase conformation classification); Kincore database, NAR 2022.
- KiSSim: Sydow et al., J. Chem. Inf. Model. 2022.
- Davis et al., Nature Biotechnology 2011 (comprehensive kinase inhibitor profiling).
- Manning et al., Science 2002 (the protein kinase complement of the human genome).
- P2Rank: Krivák and Hoksza, J. Cheminformatics 2018.
- fpocket: Le Guilloux et al., BMC Bioinformatics 2009.
- IQ-TREE 2: Minh et al., Mol. Biol. Evol. 2020.
- MMseqs2: Steinegger and Söding, Nat. Biotechnol. 2017.
- Mol\*: Sehnal et al., Nucleic Acids Research 2021.
