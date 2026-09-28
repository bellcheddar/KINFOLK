"""Phase 2c: the ESMFold2 fold queue (spec 4.4, 2.3).

One long-lived process holds one resident model. Jobs are (domain, sample
range) batches, ordered so that a partial run still yields the headline set:

  pass A  T0 controls + ABL1, 32 samples
  pass B  T1 blind + novel hold-outs, 16 samples
  pass C  one sample for every other domain (minimum viable set, gate G1)
  pass D  T2 pseudokinases + dark kinome, up to 16 samples
  pass Q6 T4 ligand-conditioned folds (apo and holo), 8 samples each
  pass E  T3 remaining domains, up to the T3 sample count

Within a pass, sequences are bucketed by 50 residues, shortest first. Each
finished batch appends to data/manifests/fold.jsonl; the annotation pool polls it.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import time
import traceback
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

from . import DATA
from .config import Manifest, log, tuning

MODELS = DATA / "esm2" / "models"
EXTRA_JOBS = DATA / "esm2" / "extra_jobs.jsonl"
STATUS = DATA / "esm2" / "queue_status.json"
SAMPLES = {"T0": 32, "T1": 16, "T2": 16, "T3": 8, "T4": 8}
PAD = 10
BATCH = 8


@dataclass
class Job:
    job_id: str
    domain_id: str
    sequence: str
    construct_start: int
    pass_name: str
    first_sample: int
    n_samples: int
    ligands: list[str] = field(default_factory=list)
    tag: str = "apo"

    @property
    def out_dir(self):
        return MODELS / self.domain_id / self.tag


def construct(row) -> tuple[str, int]:
    full = row.full_sequence
    s = max(1, row.domain_start - PAD)
    e = min(len(full), row.domain_end + PAD)
    return full[s - 1:e], s


def samples_for(tier: str) -> int:
    t = tuning()
    return int(t.get(f"samples_{tier}", SAMPLES[tier]))


def _batches(domain_id, seq, start, pass_name, lo, hi, ligands=None, tag="apo"):
    out = []
    for a in range(lo, hi, BATCH):
        n = min(BATCH, hi - a)
        out.append(Job(f"{domain_id}:{tag}:{a}-{a + n - 1}", domain_id, seq, start, pass_name, a, n,
                       ligands or [], tag))
    return out


def extra_jobs() -> list[Job]:
    if not EXTRA_JOBS.exists():
        return []
    jobs = []
    for line in EXTRA_JOBS.read_text().splitlines():
        if line.strip():
            d = json.loads(line)
            jobs += _batches(d["domain_id"], d["sequence"], d["construct_start"], d.get("pass", "Q6"),
                             0, d["n_samples"], d.get("ligands"), d.get("tag", "apo"))
    return jobs


def build_jobs() -> list[Job]:
    cat = pd.read_parquet(DATA / "catalogue" / "kinases.parquet")
    cat["construct"], cat["cstart"] = zip(*cat.apply(construct, axis=1))
    cat["bucket"] = cat.construct.str.len() // int(tuning().get("fold_bucket_width", 50))
    cat = cat.sort_values(["bucket", "domain_id"])
    jobs: list[Job] = []

    def add(rows, pass_name, lo, hi_fn):
        for _, r in rows.iterrows():
            hi = hi_fn(r)
            if hi > lo:
                jobs.extend(_batches(r.domain_id, r.construct, r.cstart, pass_name, lo, hi))

    t0, t1 = cat[cat.priority_tier == "T0"], cat[cat.priority_tier == "T1"]
    rest = cat[~cat.priority_tier.isin(["T0", "T1"])]
    add(t0, "A", 0, lambda r: samples_for("T0"))
    add(t1, "B", 0, lambda r: samples_for("T1"))
    add(rest, "C", 0, lambda r: 1)
    add(rest[rest.priority_tier == "T2"], "D", 1, lambda r: samples_for("T2"))
    q6 = [j for j in extra_jobs()]
    jobs.extend(j for j in q6 if j.pass_name != "A")
    add(rest[rest.priority_tier == "T3"], "E", 1, lambda r: samples_for("T3"))
    # Control special cases (e.g. AKT1 PH + kinase construct) go in pass A.
    jobs[0:0] = [j for j in q6 if j.pass_name == "A"]
    return jobs


def seed_for(job: Job) -> int:
    h = hashlib.sha1(f"{job.domain_id}:{job.tag}:{job.first_sample}".encode()).hexdigest()
    return int(h[:8], 16)


def estimate_seconds(length: int, n: int) -> float:
    t = tuning()
    a, b, l0 = t.get("fold_fixed_s"), t.get("fold_per_sample_s"), t.get("fold_ref_length")
    if a is None:
        return float("nan")
    return (a + b * n) * (length / l0) ** 2


def projection(jobs: list[Job] | None = None) -> dict:
    jobs = jobs if jobs is not None else build_jobs()
    done = Manifest("fold").done()
    by_pass: dict[str, float] = {}
    for j in jobs:
        if j.job_id in done:
            continue
        by_pass[j.pass_name] = by_pass.get(j.pass_name, 0) + estimate_seconds(len(j.sequence), j.n_samples)
    hours = {k: round(v / 3600, 1) for k, v in by_pass.items()}
    mvs = sum(v for k, v in by_pass.items() if k in ("A", "B", "C")) / 3600
    return {"hours_by_pass": hours, "total_hours": round(sum(by_pass.values()) / 3600, 1),
            "minimum_viable_hours": round(mvs, 1), "jobs_remaining": sum(1 for j in jobs if j.job_id not in done)}


def _write_sample(job: Job, idx: int, s) -> dict:
    job.out_dir.mkdir(parents=True, exist_ok=True)
    stem = job.out_dir / f"s{idx:02d}"
    with gzip.open(f"{stem}.cif.gz", "wt") as fh:
        fh.write(s.cif)
    np.savez_compressed(f"{stem}.npz", plddt=s.plddt.astype(np.float32),
                        pae=(s.pae.astype(np.float16) if s.pae is not None else np.zeros(0, np.float16)))
    return {"sample": idx, "mean_plddt": round(float(np.mean(s.plddt)), 2) if s.plddt.size else None,
            "ptm": round(s.ptm, 4)}


RECYCLE_JOBS = 40      # restart the worker process this often: MPS caches compiled graphs per input
                       # shape, and over ~12 h of varied lengths "other allocations" reached 156 GB
MAX_OOM_RETRIES = 2    # fresh-process OOMs before a batch is recorded as failed
EXIT_COMPLETE, EXIT_RECYCLE = 0, 3


def _is_oom(msg: str) -> bool:
    m = msg.lower()
    return "out of memory" in m or "alloc" in m or "mps backend" in m


def _release(folder) -> None:
    import gc
    gc.collect()
    if folder.device == "mps":
        folder.torch.mps.empty_cache()


def run(max_hours: float | None = None, max_jobs: int | None = None) -> int:
    """Fold until done (returns EXIT_COMPLETE) or until this worker should be recycled
    (returns EXIT_RECYCLE: after max_jobs batches, or after an out-of-memory error)."""
    from .fold import Folder

    man = Manifest("fold")
    folder = Folder()
    t = tuning()
    loops, steps = int(t.get("num_loops", 3)), int(t.get("num_sampling_steps", 100))
    log(f"fold queue: worker started on {folder.device} {folder.dtype_name}, loops {loops}, steps {steps}")
    t_start, n_jobs = time.time(), 0
    while True:
        jobs = build_jobs()
        recs = man.records()
        done = {r["id"] for r in recs if r.get("status") in ("done", "absent")}
        failed = {r["id"] for r in recs if str(r.get("status", "")).startswith("failed")}
        todo = [j for j in jobs if j.job_id not in done and j.job_id not in failed]
        if not todo:
            log("fold queue: all jobs complete")
            return EXIT_COMPLETE
        job = todo[0]
        batch = job.n_samples
        while True:
            try:
                t0 = time.time()
                samples = []
                for a in range(0, job.n_samples, batch):
                    n = min(batch, job.n_samples - a)
                    samples += folder.fold(job.sequence, samples=n, seed=seed_for(job) + a, num_loops=loops,
                                           num_sampling_steps=steps, ligands=job.ligands)
                out = [_write_sample(job, job.first_sample + i, s) for i, s in enumerate(samples)]
                man.add(job.job_id, status="done", domain_id=job.domain_id, tag=job.tag, pass_name=job.pass_name,
                        length=len(job.sequence), construct_start=job.construct_start, ligands=job.ligands,
                        seconds=round(time.time() - t0, 1), samples=out)
                break
            except Exception as e:  # noqa: BLE001
                msg = f"{type(e).__name__}: {str(e)[:160]}"
                _release(folder)
                if _is_oom(msg) and batch > 1:
                    batch = max(1, batch // 2)
                    log(f"fold queue: {job.job_id} memory error, retrying with batch {batch}")
                    continue
                if _is_oom(msg):
                    tries = sum(1 for r in recs if r["id"] == job.job_id and r.get("status") == "oom_retry") + 1
                    if tries < MAX_OOM_RETRIES:
                        man.add(job.job_id, status="oom_retry", domain_id=job.domain_id, pass_name=job.pass_name)
                        log(f"fold queue: {job.job_id} out of memory at batch 1; recycling the worker to retry")
                        _status(jobs, man)
                        return EXIT_RECYCLE
                traceback.print_exc()
                man.add(job.job_id, status=f"failed:{msg}", domain_id=job.domain_id, pass_name=job.pass_name)
                log(f"fold queue: {job.job_id} failed ({msg})")
                break
        _release(folder)
        _status(jobs, man)
        n_jobs += 1
        if max_hours and (time.time() - t_start) / 3600 > max_hours:
            log("fold queue: max_hours reached, stopping cleanly")
            return EXIT_RECYCLE
        if max_jobs and n_jobs >= max_jobs:
            log(f"fold queue: recycling worker after {n_jobs} batches")
            return EXIT_RECYCLE


def supervise() -> None:
    """Run fold workers in fresh processes until the queue is complete."""
    import subprocess
    import sys
    log(f"fold queue: supervisor started (recycle every {RECYCLE_JOBS} batches)")
    while True:
        rc = subprocess.run([sys.executable, "-c",
                             f"import sys; from kinfolk import foldqueue as q; sys.exit(q.run(max_jobs={RECYCLE_JOBS}))"]
                            ).returncode
        if rc == EXIT_COMPLETE:
            log("fold queue: supervisor finished")
            return
        if rc != EXIT_RECYCLE:
            log(f"fold queue: worker exited with code {rc}; restarting in 30 s")
            time.sleep(30)


def _status(jobs: list[Job], man: Manifest) -> None:
    recs = man.records()
    done = {r["id"] for r in recs if r["status"] == "done"}
    per_pass = {}
    for j in jobs:
        p = per_pass.setdefault(j.pass_name, [0, 0])
        p[1] += 1
        p[0] += j.job_id in done
    proj = projection(jobs)
    STATUS.write_text(json.dumps({"updated": time.strftime("%Y-%m-%d %H:%M"), "per_pass_done_total": per_pass,
                                  "failed": sum(1 for r in recs if str(r["status"]).startswith("failed")),
                                  "projection": proj}, indent=1))


if __name__ == "__main__":
    supervise()
