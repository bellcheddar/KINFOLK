"""Read config/*.toml and write the shared logs (BUILD_LOG, DECISIONS)."""
from __future__ import annotations

import datetime as dt
import json
import tomllib
from functools import lru_cache
from pathlib import Path

import tomli_w

from . import CONFIG, DATA, ROOT


@lru_cache(maxsize=1)
def tuning() -> dict:
    p = CONFIG / "tuning.toml"
    return tomllib.loads(p.read_text()) if p.exists() else {}


def update_tuning(**changes) -> None:
    t = dict(tuning())
    t.update(changes)
    (CONFIG / "tuning.toml").write_text(tomli_w.dumps(t))
    tuning.cache_clear()


def now() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M")


def log(event: str) -> None:
    with open(ROOT / "BUILD_LOG.md", "a") as fh:
        fh.write(f"- {now()}  {event}\n")


def decision(title: str, decided: str, alternatives: str, reason: str, reverse: str) -> None:
    with open(ROOT / "DECISIONS.md", "a") as fh:
        fh.write(f"## {now()}  {title}\n- Decision: {decided}\n- Alternatives: {alternatives}\n"
                 f"- Reason: {reason}\n- How to reverse: {reverse}\n\n")


class Manifest:
    """Append-only JSONL manifest for resumable stages (data/manifests/<stage>.jsonl)."""

    def __init__(self, stage: str):
        self.path = DATA / "manifests" / f"{stage}.jsonl"
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def records(self) -> list[dict]:
        if not self.path.exists():
            return []
        out = []
        for line in self.path.read_text().splitlines():
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue  # tolerate a line truncated by a kill
        return out

    def done(self) -> set[str]:
        return {r["id"] for r in self.records() if r.get("status") in ("done", "absent")}

    def add(self, item_id: str, status: str = "done", **extra) -> None:
        with open(self.path, "a") as fh:
            fh.write(json.dumps({"id": item_id, "status": status, "at": now(), **extra}) + "\n")
