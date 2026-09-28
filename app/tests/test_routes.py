"""QC: every route returns 200 (20 random kinases); every specimen JSON loads and fits the budget."""
import json
import random
import sqlite3
import sys
from pathlib import Path

import pytest

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))
DB = APP / "data" / "kinfolk.sqlite"
pytestmark = pytest.mark.skipif(not DB.exists(), reason="bundle not built")


@pytest.fixture(scope="module")
def client():
    from kinfolk_app import create_app
    return create_app().test_client()


def kinases():
    con = sqlite3.connect(DB)
    ids = [r[0] for r in con.execute("SELECT domain_id FROM kinase")]
    random.Random(13).shuffle(ids)
    return ids


@pytest.mark.parametrize("path", ["/", "/resemblance", "/rooms", "/black-sheep", "/dark", "/almanac", "/about",
                                  "/api/kinases", "/api/sites", "/api/pseudo", "/api/dark",
                                  "/static/data/tree_sequence.json"])
def test_pages(client, path):
    assert client.get(path).status_code == 200


def test_twenty_random_kinases(client):
    for k in kinases()[:20]:
        assert client.get(f"/k/{k}").status_code == 200, k
        r = client.get(f"/static/data/kinase/{k}.json")
        assert r.status_code == 200, k
        json.loads(r.data)
        assert client.get(f"/api/similar/{k}").status_code == 200
        assert client.get(f"/api/sites/{k}").status_code == 200


def test_every_specimen_json_loads_and_fits():
    big = []
    for p in (APP / "static" / "data" / "kinase").glob("*.json"):
        json.loads(p.read_text())
        if p.stat().st_size > 150_000:
            big.append((p.name, p.stat().st_size))
    assert not big, big


def test_404(client):
    assert client.get("/k/NOT_A_KINASE").status_code == 404


def test_structures_are_gzip_encoded(client):
    con = sqlite3.connect(DB)
    row = con.execute("SELECT file_path FROM structure WHERE file_path IS NOT NULL LIMIT 1").fetchone()
    if not row:
        pytest.skip("no bundled structures")
    r = client.get(f"/structures/{row[0]}")
    assert r.status_code == 200 and r.headers.get("Content-Encoding") == "gzip"
