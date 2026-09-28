"""KINFOLK field guide: Flask app serving the precomputed atlas (read-only SQLite + static JSON)."""
from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import orjson
from flask import Flask, Response, abort, g, render_template, request, send_from_directory

HERE = Path(__file__).resolve().parent
DB_PATH = Path(os.environ.get("KINFOLK_DB", HERE / "data" / "kinfolk.sqlite"))
STRUCT_DIR = HERE / "data" / "structures"

PLATES = [
    ("I", "Family Tree", "/", "tree"),
    ("II", "Specimen Card", "/k/ABL1", "specimen"),
    ("III", "Family Resemblance", "/resemblance", "resemblance"),
    ("IV", "Hidden Rooms", "/rooms", "rooms"),
    ("V", "Black Sheep", "/black-sheep", "black_sheep"),
    ("VI", "Distant Relatives", "/dark", "dark"),
    ("VII", "The Almanac", "/almanac", "almanac"),
]


def create_app() -> Flask:
    app = Flask(__name__, static_folder="static", template_folder="templates")
    app.config["JSON_SORT_KEYS"] = False
    app.jinja_env.filters["fromjson"] = lambda v: orjson.loads(v) if v else {}

    def db() -> sqlite3.Connection:
        if "db" not in g:
            g.db = sqlite3.connect(f"file:{DB_PATH}?mode=ro&immutable=1", uri=True, check_same_thread=False)
            g.db.row_factory = sqlite3.Row
        return g.db

    @app.teardown_appcontext
    def _close(_exc):
        conn = g.pop("db", None)
        if conn is not None:
            conn.close()

    def rows(sql: str, *args) -> list[dict]:
        return [dict(r) for r in db().execute(sql, args).fetchall()]

    def one(sql: str, *args) -> dict | None:
        r = db().execute(sql, args).fetchone()
        return dict(r) if r else None

    def summary(key: str, default=None):
        r = one("SELECT value_json FROM summary WHERE key = ?", key)
        return orjson.loads(r["value_json"]) if r else default

    def js(obj) -> Response:
        return Response(orjson.dumps(obj, option=orjson.OPT_SERIALIZE_NUMPY), mimetype="application/json")

    @app.context_processor
    def _ctx():
        return {"plates": PLATES, "path": request.path, "build": summary("build", {})}

    # ------------------------------------------------------------ pages
    @app.get("/")
    def tree():
        return render_template("tree.html", plate="I", title="Family Tree")

    @app.get("/k/<domain_id>")
    def specimen(domain_id: str):
        k = one("SELECT * FROM kinase WHERE domain_id = ?", domain_id)
        if k is None:
            abort(404)
        return render_template("specimen.html", plate="II", title=f"{k['gene']}", k=k)

    @app.get("/resemblance")
    def resemblance():
        return render_template("resemblance.html", plate="III", title="Family Resemblance",
                               q5=summary("q5", {}), tangle=summary("tanglegram", {}))

    @app.get("/rooms")
    def rooms():
        return render_template("rooms.html", plate="IV", title="Hidden Rooms",
                               controls=rows("SELECT * FROM control ORDER BY control_id"))

    @app.get("/black-sheep")
    def black_sheep():
        return render_template("black_sheep.html", plate="V", title="Black Sheep")

    @app.get("/dark")
    def dark():
        return render_template("dark.html", plate="VI", title="Distant Relatives")

    @app.get("/almanac")
    def almanac():
        return render_template("almanac.html", plate="VII", title="The Almanac", findings=summary("findings", {}),
                               pci=summary("pci", {}), q5=summary("q5", {}))

    @app.get("/about")
    def about():
        return render_template("about.html", plate=None, title="About", refs=summary("references", []),
                               pipeline=summary("pipeline", {}))

    # ------------------------------------------------------------ data APIs
    @app.get("/api/kinases")
    def api_kinases():
        return js(rows("SELECT * FROM kinase_tips ORDER BY domain_id"))

    @app.get("/api/similar/<domain_id>")
    def api_similar(domain_id: str):
        return js(rows("SELECT domain_b, source, kissim_sim, davis_coinhib FROM similarity "
                       "WHERE domain_a = ? ORDER BY source, kissim_sim DESC", domain_id))

    @app.get("/api/sites")
    def api_sites():
        return js(rows("SELECT s.*, k.gene, k.grp FROM site s JOIN kinase k USING(domain_id) "
                       "WHERE s.candidate = 1 OR s.known_in_pdb = 0 ORDER BY s.rank_score DESC"))

    @app.get("/api/sites/<domain_id>")
    def api_sites_domain(domain_id: str):
        return js(rows("SELECT * FROM site WHERE domain_id = ? ORDER BY recovery DESC", domain_id))

    @app.get("/api/pseudo")
    def api_pseudo():
        return js(rows("SELECT * FROM pseudo_view ORDER BY gene"))

    @app.get("/api/dark")
    def api_dark():
        return js(rows("SELECT * FROM dark_view ORDER BY gene"))

    @app.get("/api/summary/<key>")
    def api_summary(key: str):
        v = summary(key)
        if v is None:
            abort(404)
        return js(v)

    @app.get("/structures/<path:name>")
    def structures(name: str):
        """Serve pre-gzipped mmCIF with Content-Encoding so viewers receive plain text."""
        if not name.endswith(".cif.gz"):
            abort(404)
        resp = send_from_directory(STRUCT_DIR, name, mimetype="text/plain", max_age=86400)
        resp.headers["Content-Encoding"] = "gzip"
        resp.headers["Vary"] = "Accept-Encoding"
        return resp

    @app.errorhandler(404)
    def not_found(_e):
        return render_template("404.html", plate=None, title="Not found"), 404

    return app


app = create_app()
