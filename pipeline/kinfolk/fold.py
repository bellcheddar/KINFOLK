"""Resident ESMFold2 wrapper (one model per process, never reloaded per job)."""
from __future__ import annotations

import os
import time
from dataclasses import dataclass

import numpy as np

from .config import tuning

DTYPES = {"float32": "float32", "bfloat16": "bfloat16", "float16": "float16"}


def apply_env() -> None:
    for k, v in tuning().get("env", {}).items():
        os.environ.setdefault(k, str(v))
    ratio = tuning().get("mps_high_watermark_ratio")
    if ratio is not None:
        os.environ.setdefault("PYTORCH_MPS_HIGH_WATERMARK_RATIO", str(ratio))


@dataclass
class Sample:
    cif: str
    plddt: np.ndarray        # per residue, 0..100
    pae: np.ndarray | None   # residue x residue, Angstrom
    ptm: float
    seconds: float


class Folder:
    def __init__(self, device: str | None = None, dtype: str | None = None):
        apply_env()
        import torch
        from esm.models.esmfold2 import ESMFold2InputBuilder, EsmFold2Model

        t = tuning()
        self.torch = torch
        self.device = device or t.get("fold_device", "cpu")
        self.dtype_name = dtype or t.get("fold_dtype", "float32")
        self.model = EsmFold2Model.from_pretrained(
            "biohub/ESMFold2", device=self.device, dtype=getattr(torch, self.dtype_name)
        ).eval()
        self.builder = ESMFold2InputBuilder()

    def set_dtype(self, dtype: str) -> None:
        self.model = self.model.to(getattr(self.torch, dtype))
        self.dtype_name = dtype

    def fold(self, sequence: str, *, samples: int = 1, seed: int = 0, num_loops: int = 3,
             num_sampling_steps: int = 100, ligands: list[str] | None = None) -> list[Sample]:
        from esm.models.esmfold2 import LigandInput, ProteinInput, StructurePredictionInput

        chains = [ProteinInput(id="A", sequence=sequence)]
        for i, ccd in enumerate(ligands or []):
            chains.append(LigandInput(id=f"L{i}" if i else "L", ccd=[ccd]))
        spi = StructurePredictionInput(sequences=chains)
        t0 = time.perf_counter()
        with self.torch.inference_mode():
            res = self.builder.fold(self.model, spi, num_loops=num_loops,
                                    num_sampling_steps=num_sampling_steps,
                                    num_diffusion_samples=samples, seed=seed)
        if self.device == "mps":
            self.torch.mps.synchronize()
        dt = time.perf_counter() - t0
        res = res if isinstance(res, list) else [res]
        out = []
        for r in res:
            plddt = r.plddt.detach().float().cpu().numpy().reshape(-1) if r.plddt is not None else np.array([])
            if plddt.size and plddt.max() <= 1.0:
                plddt = plddt * 100.0
            pae = r.pae.detach().float().cpu().numpy().squeeze() if r.pae is not None else None
            out.append(Sample(cif=r.complex.to_mmcif(), plddt=plddt, pae=pae,
                              ptm=float(r.ptm) if r.ptm is not None else float("nan"),
                              seconds=dt / len(res)))
        return out
