# Science sanity checks

- [x] ABL1 story: PDB spatial states ['DFGin', 'DFGinter', 'DFGout', 'Unassigned']; AFDB ['DFGin']; ESMFold2 best ['DFGin']
- [x] Pseudokinase flags match Phase 1: TRIB1 DFG→SLE, ERBB3 HRD→HRN, KSR2 K→R
- [x] PCI has a CI: {"afdb": {"jsd": 0.0684, "ci95": [0.0518, 0.1017]}, "esm_single": {"jsd": 0.0785, "ci95": [0.0612, 0.1157]}, "esm_ensemble": {"jsd": 0.0803, "ci95": [0.0625, 0.1148]}}
- [x] Q5 PDB rho > 0 with CI excluding 0: rho 0.298, CI [0.217, 0.376]
- [x] Positive controls evaluated: recovered ['ABL1', 'MAP2K1', 'EGFR'] of ['ABL1', 'MAP2K1', 'PDPK1', 'EGFR', 'AKT1']
