# KINFOLK findings

Headline results for the six questions in the build spec. Intervals are 95% bootstrap intervals over kinases (1,000 resamples) unless stated. Figures are live in the app at the paths given.

## Q1: Motif fidelity

Against the nearest experimental state, AlphaFold DB motifs deviate by a median of g.l 1.02 Å, aC 0.57 Å, xDFG 0.60 Å, a.l 0.70 Å; ESMFold2's best sample by g.l 1.11 Å, aC 0.63 Å, xDFG 0.68 Å, a.l 0.89 Å. The least faithful motif is the g.l region, as expected for the most mobile elements.

Figure: `/almanac#q1`

## Q2: Conformational prior (Prior Collapse Index)

Across 289 kinases, the Prior Collapse Index (Jensen–Shannon divergence from the PDB state distribution) is 0.068 for AlphaFold DB (95% CI 0.052 to 0.102), 0.079 for the best ESMFold2 sample (0.061 to 0.116) and 0.080 for the ESMFold2 ensemble (0.062 to 0.115). The ensemble does not recover more of the PDB's state diversity than a single model.

Figure: `/almanac#q2 (PCI slider)`

## Q3: Novel sequences and structures

ESMFold2's median pocket RMSD is 1.21 Å on training-era structures and 0.95 Å on blind, sequence-novel ones (n 10); AlphaFold DB's is 0.96 and 0.96 Å respectively; on the dark kinome the two predictors' pockets differ by a median 0.53 Å; AlphaFold DB restores the degraded catalytic residue to canonical geometry in 30 of 81 pseudokinase motifs; ESMFold2 restores the degraded catalytic residue to canonical geometry in 36 of 85 pseudokinase motifs.

Figure: `/almanac#q3, /dark, /black-sheep`

## Q4: Cryptic pockets

Apo ESMFold2 ensembles recover 3 of 5 positive-control allosteric sites (ABL1, MAP2K1, EGFR); the single AlphaFold DB model recovers 2. Kinome-wide, 24 candidate cryptic sites pass the filters across 18 kinases.

Figure: `/rooms, /almanac#q4`

## Q5: Selectivity

Over 175 kinases in the Davis panel, KiSSim-lite pocket similarity correlates with inhibitor co-inhibition at Spearman ρ 0.30 for PDB pockets (95% CI 0.22 to 0.38). AlphaFold DB pockets reach 1.04 of the PDB correlation (95% CI 0.89 to 1.22). ESMFold2 pockets reach 1.04 of the PDB correlation (95% CI 0.89 to 1.23).

Figure: `/resemblance, /almanac#q5`

## Q6: Induced fit

On 30 post-cutoff type II complexes, folding with the ligand raises the share of DFG-out samples from 0% (apo) to 100% (holo) (paired Wilcoxon p = 4.3e-08); the median ligand heavy-atom RMSD to the crystal is 1.0 Å over the 15 complexes where ligand atoms could be matched by name.

Figure: `/almanac#q6`

### Q3 detail: pocket Cα RMSD to the best-resolution chain in each set (median, IQR, n)

| Set | AlphaFold DB | ESMFold2 best sample |
|---|---|---|
| Training era (before 2018-04-30) | 0.96 (0.61 to 1.83), n 258 | 1.21 (0.78 to 1.91), n 259 |
| Blind for ESMFold2 (after 2021-09-30) | 1.12 (0.72 to 1.78), n 155 | 1.20 (0.87 to 1.85), n 156 |
| Blind and sequence-novel | 0.96 (0.66 to 2.23), n 10 | 0.95 (0.67 to 1.91), n 10 |

## Caveats

- Pocket fingerprints are KiSSim-lite, a reimplementation following the KiSSim design; kissim itself could not be installed.
- ESMFold2 ran with 3 recycling loops and 100 sampling steps in float32 on Apple MPS (see DECISIONS.md).
- ABL1 with imatinib has no post-cutoff entry, so it is reported separately as "memorisation possible".
- PLK5, TAF1_C and TAF1_N could not be placed in KLIFS numbering and are excluded from pocket metrics.
- Neighbour joining needs a complete distance matrix, so the pocket trees also drop domains whose fingerprints share fewer than 50 KLIFS rows with the rest: EIF2AK4_C from the PDB tree, and PLK5, STK40, TAF1_C and TAF1_N from the ESMFold2 tree (see data/phylo/pocket_tree_info.json).
