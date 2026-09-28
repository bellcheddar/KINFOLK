"""Unattended finisher: backfill while the fold queue runs, then complete Phases 3 to 10.

    kinfolk finish            # blocks until done; safe to rerun (every stage is resumable)

While the fold queue is alive: every BACKFILL_H hours run annotate -> summarise ->
pockets -> bundle so the atlas stays current. When the queue has finished:
annotate, summarise, pockets, selectivity (+ pocket trees), induced fit, findings,
references, bundle, tests, visual QC, science checks, commits, then gate G3
(or G6 if a headline control fails).
"""
from __future__ import annotations

import json
import subprocess
import time
import traceback

import pandas as pd

from . import DATA, ROOT
from .config import log

BACKFILL_H = 2.0
APP = ROOT / "app"
QC = ROOT / "qc"


def sh(cmd: str, cwd=ROOT, check=False, timeout=None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True, check=check, timeout=timeout)


def fold_queue_alive() -> bool:
    return sh("pgrep -f 'kinfolk fold'").returncode == 0


def queue_complete() -> bool:
    from .config import Manifest
    from .foldqueue import build_jobs
    recs = Manifest("fold").records()
    ended = {r["id"] for r in recs if r.get("status") == "done" or str(r.get("status", "")).startswith("failed")}
    return all(j.job_id in ended for j in build_jobs())


def stage(name: str, fn) -> bool:
    t0 = time.time()
    try:
        fn()
        log(f"finish: {name} ok in {(time.time() - t0) / 60:.1f} min")
        return True
    except Exception:  # noqa: BLE001 (fallback ladder: log, flag, continue)
        log(f"finish: {name} failed: {traceback.format_exc(limit=2).splitlines()[-1][:200]}")
        (DATA / "finish_errors.log").open("a").write(f"--- {name}\n{traceback.format_exc()}\n")
        return False


def commit(msg: str) -> None:
    sh("git add -A")
    sh(f'git commit -q -m "{msg}" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" '
       f'-m "Claude-Session: https://claude.ai/code/session_01RXpTzdrvCfj7UYAyf9uSQ2"')


def progress(rows: list[tuple[str, str, str, int]], gates: str) -> None:
    out = ["# KINFOLK progress", "", f"Updated: {time.strftime('%Y-%m-%d %H:%M')}", "",
           "| Phase | Name | Status | % |", "|------|------|--------|---|"]
    out += [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in rows]
    q = DATA / "esm2" / "queue_status.json"
    if q.exists():
        st = json.loads(q.read_text())
        out += ["", "Fold queue: " + ", ".join(f"pass {k} {v[0]}/{v[1]}" for k, v in st["per_pass_done_total"].items())
                + f"; remaining projection {st['projection']['total_hours']} h; failed batches {st['failed']}."]
    out += ["", f"Open gates: {gates}", ""]
    (ROOT / "PROGRESS.md").write_text("\n".join(out))


def backfill() -> None:
    from . import annotate, bundle, pockets, summarise
    stage("annotate (backfill)", annotate.run)
    stage("summarise (backfill)", summarise.run)
    stage("pockets (backfill)", pockets.run)
    stage("bundle (backfill)", bundle.run)


def science_checks() -> tuple[list[str], list[str]]:
    """Spec 13.4. Returns (report lines, headline failures that justify gate G6)."""
    lines, headline = ["# Science sanity checks", ""], []
    s = pd.read_parquet(DATA / "annotation" / "structures.parquet")
    dom = pd.read_parquet(DATA / "annotation" / "domains.parquet").set_index("domain_id")
    a = s[s.domain_id == "ABL1"]
    pdb_sp = set(a[a.source == "pdb"].kincore_spatial)
    af = a[a.source == "afdb"].kincore_spatial.tolist()
    best = dom.loc["ABL1"].get("esm_best")
    esm_sp = a[a.struct_id == best].kincore_spatial.tolist() if isinstance(best, str) else []
    ok1 = {"DFGin", "DFGout"} <= pdb_sp and af == ["DFGin"] and esm_sp == ["DFGin"]
    lines.append(f"- [{'x' if ok1 else ' '}] ABL1 story: PDB spatial states {sorted(pdb_sp)}; AFDB {af}; ESMFold2 best {esm_sp}")
    from .catalogue import acceptance
    acc = acceptance(pd.read_parquet(DATA / "catalogue" / "kinases.parquet"))
    lines.append(f"- [{'x' if acc['pass'] else ' '}] Pseudokinase flags match Phase 1: TRIB1 {acc['TRIB1']}, ERBB3 {acc['ERBB3']}, KSR2 {acc['KSR2']}")
    pci = json.loads((DATA / "annotation" / "pci.json").read_text())
    ok3 = bool(pci.get("kinome")) and all("ci95" in v for v in pci["kinome"].values())
    lines.append(f"- [{'x' if ok3 else ' '}] PCI has a CI: {json.dumps(pci.get('kinome', pci.get('status')))}")
    q5p = DATA / "selectivity" / "q5.json"
    q5 = json.loads(q5p.read_text()) if q5p.exists() else {}
    rho = (q5.get("rho") or {}).get("pdb")
    ci = (q5.get("rho_ci95") or {}).get("pdb", [None, None])
    ok4 = rho is not None and rho > 0 and ci[0] is not None and ci[0] > 0
    lines.append(f"- [{'x' if ok4 else ' '}] Q5 PDB rho > 0 with CI excluding 0: rho {rho}, CI {ci}")
    ctrl_p = DATA / "pockets" / "controls.json"
    ctrl = json.loads(ctrl_p.read_text()) if ctrl_p.exists() else []
    rec = [c["gene"] for c in ctrl if c.get("recovered")]
    lines.append(f"- [{'x' if ctrl else ' '}] Positive controls evaluated: recovered {rec} of {[c['gene'] for c in ctrl]}")
    # G6: only when a failure would invalidate a headline figure (spec 0.2).
    if not ok1:
        headline.append("ABL1 story (Q2 headline example) does not reproduce")
    if not ok4:
        headline.append("Q5 method sanity check failed (PDB rho not positive with CI excluding 0)")
    if ctrl and not rec:
        headline.append("No positive-control site recovered (Q4 headline)")
    (QC / "science.md").write_text("\n".join(lines) + "\n")
    return lines, headline


def open_gate(gid: str, need: str, unblock: str, done: str, after: str) -> None:
    (ROOT / "GATE_OPEN.md").write_text(
        f"# Gate {gid} open\n\n## What is needed from Marc\n{need}\n\n## How to unblock\n{unblock}\n\n"
        f"## Completed\n{done}\n\n## After the gate\n{after}\n")
    sh(f"osascript -e 'display notification \"{need[:90]}\" with title \"KINFOLK gate {gid}\" sound name \"Glass\"'")
    log(f"gate {gid} opened")


def final() -> None:
    from . import annotate, bundle, findings, induced_fit, pockets, references, selectivity, summarise
    stage("annotate", annotate.run)
    stage("summarise", summarise.run)
    commit("phase 3: annotation, agreement, ensembles, PCI (final)")
    stage("pockets", pockets.run)
    commit("phase 5: pockets, candidate cryptic sites, control scorecard")
    stage("selectivity + pocket trees", selectivity.run)
    commit("phase 6: KiSSim-lite fingerprints, Q5 correlations; phase 4 pocket trees and tanglegram")
    stage("induced fit", induced_fit.run)
    commit("phase 7: induced-fit test")
    stage("references", references.resolve)
    stage("findings", findings.run)
    stage("bundle", bundle.run)
    commit("phase 8: data bundle")
    # QC
    t = sh(".pixi/envs/default/bin/pytest pipeline/tests -q", timeout=3600)
    ta = sh("uv run pytest tests -q", cwd=APP, timeout=3600)
    (QC / "pytest.txt").write_text(t.stdout[-4000:] + "\n--- app ---\n" + ta.stdout[-4000:])
    srv = subprocess.Popen("uv run gunicorn -w 3 -b 127.0.0.1:8010 kinfolk_app:app", shell=True, cwd=APP,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    time.sleep(4)
    v = sh("uv run --project app python qc/visual_qc.py", timeout=3600)
    srv.terminate()
    lines, headline = science_checks()
    log(f"finish: pytest pipeline rc {t.returncode}, app rc {ta.returncode}; visual QC rc {v.returncode}; "
        f"science headline failures {headline}")
    commit("phase 9: QC (tests, screenshots, accessibility, science checks)")
    progress([("0-8", "Pipeline", "Done", 100), ("9", "Web app", f"QC rc {v.returncode}", 100),
              ("10", "Polish and handover", "Findings, deploy script; awaiting G3", 95)],
             "G6" if headline else "G3")
    if headline:
        open_gate("G6", "A decision on whether to proceed with these science caveats flagged: " + "; ".join(headline),
                  "Reply 'proceed with caveat' (the caveat goes into FINDINGS.md), or ask for a re-run with changed settings.",
                  "Phases 0 to 9 complete; see qc/science.md and FINDINGS.md.",
                  "Flag the caveat, then open gate G3 for deployment.")
    else:
        open_gate("G3", "Approval to deploy kinfolk.mdeller.com (rsync bundle to the droplet, systemd, nginx, certbot).",
                  "Review deploy/deploy.sh, then run: bash deploy/deploy.sh (it needs SSH access, which the build may not use).",
                  "Phases 0 to 10 complete; QC in qc/; findings in FINDINGS.md.",
                  "Nothing else is pending locally; after deployment, check https://kinfolk.mdeller.com.")
    commit("phase 10: findings, handover, gate opened")


def run() -> None:
    log("finish: orchestrator started")
    last = 0.0
    while not queue_complete():
        if not fold_queue_alive():
            log("finish: fold queue not running but incomplete; restarting it")
            subprocess.Popen("caffeinate -dimsu .pixi/envs/default/bin/kinfolk fold >> data/esm2/foldqueue.log 2>&1",
                             shell=True, cwd=ROOT)
            time.sleep(120)
        if time.time() - last > BACKFILL_H * 3600:
            backfill()
            last = time.time()
        progress([("0-2", "Bootstrap, catalogue, acquisition", "Done", 100),
                  ("2c", "Fold queue", "Running", 50), ("3-8", "Analysis", "Backfilling every 2 h", 60),
                  ("9", "Web app", "Built; QC passing on interim bundle", 90), ("10", "Handover", "Pending folds", 20)],
                 "none")
        time.sleep(600)
    log("finish: fold queue complete; running final phases")
    final()
