"""kinfolk CLI: `kinfolk run --phase N`, `kinfolk status`, `kinfolk fold`."""
from __future__ import annotations

import argparse
import json

from . import DATA, ROOT
from .config import Manifest


def phase(n: int) -> None:
    if n == 1:
        from . import catalogue
        df = catalogue.build()
        print(catalogue.acceptance(df))
    elif n == 2:
        from . import acquire
        st = acquire.acquire_pdb()
        acquire.assign_t1(st)
        acquire.acquire_afdb()
    elif n == 3:
        from . import annotate, summarise
        annotate.run()
        summarise.run()
    elif n == 4:
        from . import phylo
        phylo.sequence_tree()
    elif n == 5:
        from . import pockets
        pockets.run()
    elif n == 6:
        from . import selectivity
        selectivity.run()
    elif n == 7:
        from . import induced_fit
        induced_fit.run()
    elif n == 8:
        from . import bundle
        bundle.run()
    else:
        raise SystemExit(f"phase {n} is not a pipeline phase (0: tools/, 9-10: app/ and deploy/)")


def status() -> None:
    print((ROOT / "PROGRESS.md").read_text() if (ROOT / "PROGRESS.md").exists() else "no PROGRESS.md")
    q = DATA / "esm2" / "queue_status.json"
    if q.exists():
        print(json.dumps(json.loads(q.read_text()), indent=1))
    for m in sorted((DATA / "manifests").glob("*.jsonl")):
        recs = Manifest(m.stem).records()
        failed = sum(1 for r in recs if str(r.get("status", "")).startswith("failed"))
        print(f"{m.stem:24s} {len(recs):6d} records, {failed} failed")


def main() -> None:
    ap = argparse.ArgumentParser(prog="kinfolk")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--phase", type=int, required=True)
    sub.add_parser("status")
    sub.add_parser("finish")
    f = sub.add_parser("fold")
    f.add_argument("--max-hours", type=float, default=None)
    a = ap.parse_args()
    if a.cmd == "run":
        phase(a.phase)
    elif a.cmd == "status":
        status()
    elif a.cmd == "finish":
        from . import finish
        finish.run()
    elif a.cmd == "fold":
        from . import foldqueue
        if a.max_hours:
            foldqueue.run(a.max_hours)
        else:
            foldqueue.supervise()


if __name__ == "__main__":
    main()
