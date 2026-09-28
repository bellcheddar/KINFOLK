"""Phase 0 ESMFold2 smoke test on ABL1 (spec 2.2 fold_dtype, 4.4 budgeting).

1. Fold the ABL1 kinase domain (+/-10 residue padding) in float32 as reference.
2. Try bfloat16 then float16; keep the fastest dtype whose pocket-core Ca RMSD
   to the float32 reference is < 0.5 A.
3. Compare 50 vs 100 sampling steps in the chosen dtype.
4. Time an 8-sample ensemble to measure the per-sample marginal cost.
Falls back from MPS to CPU if MPS fails. Results: data/esm2/smoke/smoke.json.
"""
from __future__ import annotations

import json
import sys
import time
import traceback

import requests

from kinfolk import DATA
from kinfolk import klifs
from kinfolk.config import decision, log, tuning, update_tuning
from kinfolk.structure import atom_table, read_structure, superposed_rmsd

OUT = DATA / "esm2" / "smoke"
LOOPS, STEPS = 3, 100
RMSD_TOL = 0.5
KNOWN_BAD = {
    "mps": {
        "bfloat16": "process abort in MPSNDArrayMatrixMultiplication (destination and accumulator dtypes differ); "
                    "observed 2026-09-27 20:10, esm 3.4.1.post1, torch 2.11.0",
        "float16": "not supported: ESMC precision options are bf16 and fp32 only, and trunk and ESMC must match",
    },
}


def abl1_construct() -> tuple[str, str, int]:
    seq = "".join(requests.get("https://rest.uniprot.org/uniprotkb/P00519.fasta", timeout=60)
                  .text.split("\n")[1:])
    pocket = klifs.get("kinase_information", kinase_ID=392)[0]["pocket"]
    start, end = max(1, 242 - 10), min(len(seq), 493 + 10)
    return seq[start - 1:end], pocket, start


def core_rmsd(cif_a: str, cif_b: str, mapping: dict[int, int]) -> float:
    a, b = atom_table(read_structure(cif_a)), atom_table(read_structure(cif_b))
    keep = [p for p in klifs.CORE_POSITIONS if p in mapping and mapping[p] in a and mapping[p] in b]
    import numpy as np
    return superposed_rmsd(np.array([a[mapping[p]] for p in keep]), np.array([b[mapping[p]] for p in keep]))


def run(device: str) -> dict:
    from kinfolk.fold import Folder

    seq, pocket, start = abl1_construct()
    mapping = klifs.map_pocket_to_sequence(pocket, seq)  # residue numbers 1..len(construct)
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {"device": device, "length": len(seq), "construct_start": start,
                    "num_loops": LOOPS, "num_sampling_steps": STEPS, "dtypes": {}}

    t0 = time.perf_counter()
    f = Folder(device=device, dtype="float32")
    report["load_seconds"] = round(time.perf_counter() - t0, 1)
    log(f"smoke: ESMFold2 loaded on {device} float32 in {report['load_seconds']} s")
    f.fold(seq[:60], num_sampling_steps=10)  # warm-up (kernel compilation, caches)

    ref = f.fold(seq, seed=0, num_loops=LOOPS, num_sampling_steps=STEPS)[0]
    (OUT / "abl1_float32.cif").write_text(ref.cif)
    report["dtypes"]["float32"] = {"seconds": round(ref.seconds, 1), "core_rmsd": 0.0,
                                   "mean_plddt": round(float(ref.plddt.mean()), 1), "ptm": round(ref.ptm, 3)}
    log(f"smoke: float32 ABL1 fold {ref.seconds:.1f} s, pLDDT {ref.plddt.mean():.1f}")

    # Reduced-precision trials that cannot run in-process are recorded, not retried:
    # an MPS assertion aborts the interpreter, so it cannot be caught.
    for dt, why in KNOWN_BAD.get(device, {}).items():
        report["dtypes"][dt] = {"error": why}
    if device == "mps":
        for dt in [d for d in ("bfloat16", "float16") if d not in KNOWN_BAD.get(device, {})]:
            try:
                f.set_dtype(dt)
                f.fold(seq[:60], num_sampling_steps=10)
                s = f.fold(seq, seed=0, num_loops=LOOPS, num_sampling_steps=STEPS)[0]
                (OUT / f"abl1_{dt}.cif").write_text(s.cif)
                r = core_rmsd(s.cif, ref.cif, mapping)
                report["dtypes"][dt] = {"seconds": round(s.seconds, 1), "core_rmsd": round(r, 3),
                                        "mean_plddt": round(float(s.plddt.mean()), 1), "ptm": round(s.ptm, 3)}
                log(f"smoke: {dt} ABL1 fold {s.seconds:.1f} s, core RMSD to float32 {r:.2f} A")
            except Exception as e:  # noqa: BLE001
                report["dtypes"][dt] = {"error": f"{type(e).__name__}: {e}"[:300]}
                log(f"smoke: {dt} failed ({type(e).__name__})")

    ok = {k: v for k, v in report["dtypes"].items() if "seconds" in v and v["core_rmsd"] < RMSD_TOL}
    best = min(ok, key=lambda k: ok[k]["seconds"])
    report["chosen_dtype"] = best
    f.set_dtype(best)
    base = f.fold(seq, seed=0, num_loops=LOOPS, num_sampling_steps=STEPS)[0]

    s50 = f.fold(seq, seed=0, num_loops=LOOPS, num_sampling_steps=50)[0]
    report["steps50"] = {"seconds": round(s50.seconds, 1),
                         "core_rmsd_vs_100": round(core_rmsd(s50.cif, base.cif, mapping), 3)}

    ens = f.fold(seq, seed=1, num_loops=LOOPS, num_sampling_steps=STEPS, samples=8)
    report["ensemble8_seconds_total"] = round(sum(s.seconds for s in ens), 1)
    report["single_seconds"] = round(base.seconds, 1)
    for i, s in enumerate(ens):
        (OUT / f"abl1_ens_{i}.cif").write_text(s.cif)
    return report


def main() -> int:
    tuning()
    report = None
    for device in ("mps", "cpu"):
        try:
            report = run(device)
            break
        except Exception as e:  # noqa: BLE001
            log(f"smoke: device {device} failed: {type(e).__name__}: {str(e)[:200]}")
            traceback.print_exc()
    if report is None:
        print("SMOKE FAILED ON MPS AND CPU")
        return 2
    (OUT / "smoke.json").write_text(json.dumps(report, indent=2))
    per_sample = max(0.0, (report["ensemble8_seconds_total"] - report["single_seconds"]) / 7)
    fixed = max(0.0, report["single_seconds"] - per_sample)
    update_tuning(fold_device=report["device"], fold_dtype=report["chosen_dtype"],
                  fold_dtype_status="set by smoke test", num_loops=LOOPS, num_sampling_steps=STEPS,
                  fold_fixed_s=round(fixed, 2), fold_per_sample_s=round(per_sample, 2),
                  fold_ref_length=report["length"])
    decision("Fold precision", f"{report['chosen_dtype']} on {report['device']}",
             "float32, bfloat16, float16",
             f"Fastest dtype reproducing the float32 ABL1 pocket core to < {RMSD_TOL} A: "
             + json.dumps(report["dtypes"]),
             "Edit fold_dtype in config/tuning.toml")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
