<p align="center"><img src="kinfolk_icon_300.png" width="160" alt="KINFOLK icon: an open field guide with a kinase family tree"></p>

# 🧬 KINFOLK

> **A kinome-wide field guide to how ESMFold2, AlphaFold DB and the PDB see every human kinase pocket.**

![python](https://img.shields.io/badge/python-3.12-3776AB?logo=python&logoColor=white) ![flask](https://img.shields.io/badge/flask-3.1-467FF7?logo=flask&logoColor=white) ![ESMFold2](https://img.shields.io/badge/ESMFold2-7B%20on%20Apple%20MPS-467FF7) ![domains](https://img.shields.io/badge/kinase%20domains-522-00897B) ![status](https://img.shields.io/badge/status-pre--deployment-fcb900) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller%2C%20D.Phil.-1C244B)

<table>
<tr>
<td>🌐 <b>Website</b></td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️ <b>Contact</b></td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙 <b>GitHub</b></td><td><a href="https://github.com/bellcheddar" target="_blank" rel="noopener noreferrer">@bellcheddar</a></td>
</tr>
</table>

---

**K**inase **IN**ventory of **F**olds, **O**ccluded **L**igand sites and **K**ryptic pockets. KINFOLK compares three views of every human kinase catalytic domain: experimental structures from the PDB (via KLIFS), the single AlphaFold DB model, and an ensemble of ESMFold2 diffusion samples. It measures KLIFS pocket motifs, conformational states, novel and dark kinases, cryptic pockets, selectivity and induced fit, and serves the results as a precomputed, interactive field guide at `kinfolk.mdeller.com`.

Why it matters: structure predictors are now the default starting point for kinase drug design, yet a single predicted model tends to show one conformation, and it is rarely clear how far to trust a predicted pocket for a kinase with no crystal structure. KINFOLK quantifies where predicted pockets agree with experiment, where they collapse the kinome onto one active-like state, and where a diffusion ensemble recovers the missing states or opens hidden allosteric rooms. It is useful for: choosing which predicted model to dock into, triaging dark and pseudokinases, finding candidate allosteric sites, and anticipating off-target kinases from pocket similarity.

## ✨ Features

| Plate | Route | What it shows |
|---|---|---|
| I. Family Tree | `/` | Phylogram, cladogram or radial tree of 522 catalytic domains; three agreement marks per tip (PDB ring, AlphaFold DB square, ESMFold2 triangle) |
| II. Specimen Card | `/k/<domain_id>` | Observer's note, field marks, conformational compass (Kincore D1/D2), 85-position KLIFS ribbon, three-way Mol* overlay with two-way linking |
| III. Family Resemblance | `/resemblance` | Sequence tree against pocket tree (tanglegram), off-target finder, selectivity correlations |
| IV. Hidden Rooms | `/rooms` | Positive-control scorecard, candidate cryptic sites, site viewer, subpocket sunburst |
| V. Black Sheep | `/black-sheep` | Pseudokinases: does a predictor "repair" the degraded catalytic motif? |
| VI. Distant Relatives | `/dark` | The dark kinome, with predictor agreement as a confidence signal |
| VII. The Almanac | `/almanac` | One section per question (Q1 to Q6), including the Prior Collapse Index slider |

Every view keeps its state in the URL (for example `/?k=ABL1&tree=phylo&tips=agree`), so any view can be shared.

## 🎯 Questions

| | Question | Headline metric |
|---|---|---|
| Q1 | How faithful is each KLIFS motif? | Cα RMSD per motif to the nearest PDB state, pocket-core frame |
| Q2 | Do predictors collapse the kinome onto BLAminus? | Prior Collapse Index: Jensen–Shannon divergence from the PDB state distribution, bootstrap CI |
| Q3 | How do models behave on dark, pseudo and post-cutoff kinases? | Blind and blind + novel pocket RMSD; motif-repair calls |
| Q4 | Do apo ESMFold2 ensembles find known allosteric sites? | Positive-control recovery (Jaccard ≥ 0.4); candidate cryptic sites |
| Q5 | Do predicted pockets predict selectivity? | Spearman ρ against Davis 2011 co-inhibition, as a fraction of the PDB value |
| Q6 | Does a type II ligand pull DFG out? | DFG-out fraction, apo against holo ESMFold2 folds |

Results, with numbers and confidence intervals, are in [`FINDINGS.md`](FINDINGS.md).

## 📋 Requirements

| Component | Version | Notes |
|---|---|---|
| macOS on Apple Silicon | 15 or later | Built on a Mac Studio (M2 Ultra, 128 GB) |
| pixi | 0.6x | Compute environment (`pixi.toml`, conda-forge and bioconda) |
| uv | 0.9x | App environment (`app/pyproject.toml`) |
| ESMFold2 weights | `biohub/ESMFold2` | 26.7 GB, fetched from Hugging Face on first run |
| Google Chrome | any recent | Used by the Playwright QC (never `playwright install`) |

## 🚀 Usage

```bash
# compute environment
pixi install && bash tools/post_install.sh && bash tools/fetch_vendor.sh
pixi run probe                      # hardware probe -> config/hardware.toml, tuning.toml
pixi run smoke                      # ESMFold2 smoke test on ABL1, fold-time projection

# pipeline (every stage is resumable; manifests in data/manifests/)
pixi run kinfolk run --phase 1      # catalogue
pixi run kinfolk run --phase 2      # PDB via KLIFS, AlphaFold DB, blind and novel flags
pixi run kinfolk fold               # ESMFold2 fold queue (long-running)
pixi run kinfolk finish             # backfill while folding, then phases 3 to 10
pixi run kinfolk status

# app
cd app && uv sync && uv run gunicorn -w 2 -b 127.0.0.1:8010 kinfolk_app:app
```

| Command | Purpose |
|---|---|
| `kinfolk run --phase N` | Run one pipeline phase (1 to 8) |
| `kinfolk fold [--max-hours H]` | Start or resume the fold queue |
| `kinfolk finish` | Unattended backfill and completion, ending at deployment gate G3 |
| `kinfolk status` | Progress, fold-queue projection and manifest counts |

## 📊 Output

| Path | Contents |
|---|---|
| `app/data/kinfolk.sqlite` | Read-only atlas (kinase, structure, motif_dev, ribbon, ligand, site, similarity, tree, control, induced_fit, summary) |
| `app/static/data/kinase/*.json` | One specimen file per domain (each under 150 KB) |
| `app/data/structures/` | Domain-cropped mmCIF (gzip) in each domain's pocket-core frame |
| `FINDINGS.md` | Headline results for Q1 to Q6 |
| `DECISIONS.md`, `BUILD_LOG.md`, `PROGRESS.md` | Every judgement call, event and phase status |
| `qc/` | Tests, screenshots, accessibility and science checks |
| `deploy/` | `deploy.sh`, systemd unit, nginx server block (not run by the build) |

## 🔧 Tech stack

| Layer | Tools |
|---|---|
| Structure prediction | ESMFold2 (Biohub `esm` 3.4.1), float32 on Apple MPS |
| Structures and annotation | KLIFS, RCSB, AlphaFold DB, gemmi, biotite, Kincore rules (ported) |
| Sequence and phylogeny | HMMER (Pfam PF00069), MMseqs2, IQ-TREE 3, scikit-bio |
| Pockets | P2Rank 2.5.1, fpocket 4.2 |
| Selectivity | KiSSim-lite fingerprints, Davis 2011 (via PyTDC) |
| App | Flask, gunicorn, SQLite, D3 7.9, Mol* 5.11 (vendored, checksummed) |

## 📄 Licence

Code: MIT. Derived data follow the terms of their sources (PDB: CC0; AlphaFold DB: CC BY 4.0; KLIFS and Kincore: as published; ESMFold2 weights: MIT). Software and reference DOIs are listed on the app's About page.

---

## 👤 Author

**Marc C. Deller, D.Phil.**  
Structural biologist & drug discovery scientist  

<table>
<tr>
<td>🌐</td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️</td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙</td><td><a href="https://github.com/bellcheddar" target="_blank" rel="noopener noreferrer">github.com/bellcheddar</a></td>
</tr>
</table>
