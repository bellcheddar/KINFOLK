"""Phase 0 hardware probe (spec section 2.1 and 2.2).

Writes config/hardware.toml and config/tuning.toml. Every downstream
parallelism setting is read from tuning.toml; nothing is hard-coded.
fold_device and fold_dtype are provisional here and are finalised by
tools/smoke_fold.py once the ABL1 smoke test has run.
"""
from __future__ import annotations

import datetime as dt
import re
import shutil
import subprocess
import tomllib
from pathlib import Path

import tomli_w

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"


def sh(cmd: list[str]) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def gpu_cores() -> int:
    out = sh(["system_profiler", "SPDisplaysDataType"])
    m = re.search(r"Total Number of Cores:\s*(\d+)", out)
    return int(m.group(1)) if m else 0


def mps_available() -> bool:
    try:
        import torch

        return bool(torch.backends.mps.is_available())
    except Exception:
        return False


def probe() -> dict:
    usage = shutil.disk_usage(ROOT)
    return {
        "probed_at": dt.datetime.now().isoformat(timespec="seconds"),
        "chip": sh(["sysctl", "-n", "machdep.cpu.brand_string"]),
        "performance_cores": int(sh(["sysctl", "-n", "hw.perflevel0.physicalcpu"]) or 0),
        "efficiency_cores": int(sh(["sysctl", "-n", "hw.perflevel1.physicalcpu"]) or 0),
        "unified_memory_gb": round(int(sh(["sysctl", "-n", "hw.memsize"]) or 0) / 2**30, 1),
        "gpu_cores": gpu_cores(),
        "free_disk_gb": round(usage.free / 2**30, 1),
        "mps_available": mps_available(),
    }


def derive(hw: dict, previous: dict | None) -> dict:
    p = hw["performance_cores"]
    cpu_workers = max(2, p - 2)
    tuning = {
        "cpu_workers": cpu_workers,
        "io_workers": 16,
        "fold_device": "mps" if hw["mps_available"] else "cpu",
        "fold_dtype": "bfloat16",
        "fold_dtype_status": "provisional (set by smoke test)",
        "fold_instances": 1,
        "fold_bucket_width": 50,
        "iqtree_threads": f"-T AUTO --threads-max {p}",
        "p2rank_threads": cpu_workers,
        "memory_headroom_gb": 24,
        "max_compute_hours": 96,
        "env": {
            "PYTORCH_ENABLE_MPS_FALLBACK": "1",
            "TOKENIZERS_PARALLELISM": "false",
            "OMP_NUM_THREADS": str(cpu_workers),
        },
    }
    # Keep decisions made by the smoke test across re-probes.
    if previous:
        for key in ("fold_device", "fold_dtype", "fold_dtype_status", "num_sampling_steps",
                    "num_loops", "max_compute_hours", "mps_high_watermark_ratio"):
            if key in previous:
                tuning[key] = previous[key]
    return tuning


def main() -> None:
    CONFIG.mkdir(exist_ok=True)
    hw = probe()
    tuning_path = CONFIG / "tuning.toml"
    previous = tomllib.loads(tuning_path.read_text()) if tuning_path.exists() else None
    tuning = derive(hw, previous)
    (CONFIG / "hardware.toml").write_text(tomli_w.dumps(hw))
    tuning_path.write_text(tomli_w.dumps(tuning))
    print(tomli_w.dumps(hw))
    print(tomli_w.dumps(tuning))


if __name__ == "__main__":
    main()
