"""Visual, accessibility and timing QC (spec 13.2, 13.3, 13.5).

Uses the system Google Chrome through Playwright (channel="chrome"); never runs `playwright install`.
Usage: uv run --project app python qc/visual_qc.py [base_url]
Writes qc/screens/*.png and qc/visual_qc.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

try:
    from axe_playwright_python.sync_playwright import Axe
except ImportError:  # pragma: no cover
    Axe = None

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8010"
OUT = Path(__file__).resolve().parent
PLATES = [("I", "/"), ("II", "/k/ABL1"), ("III", "/resemblance"), ("IV", "/rooms"), ("V", "/black-sheep"),
          ("VI", "/dark"), ("VII", "/almanac"), ("About", "/about")]
VIEWPORTS = {"desktop": (1440, 900), "mobile": (390, 844)}


def check_page(page, plate: str, path: str, vp: str) -> dict:
    errors: list[str] = []
    page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
    page.on("pageerror", lambda e: errors.append(str(e)))
    resp = page.goto(BASE + path, wait_until="networkidle", timeout=60000)
    page.wait_for_timeout(2500 if plate in ("II", "IV") else 1200)
    res = {"plate": plate, "path": path, "viewport": vp, "status": resp.status if resp else None}
    res["horizontal_scroll"] = page.evaluate("document.documentElement.scrollWidth > window.innerWidth + 1")
    res["fcp_ms"] = page.evaluate(
        "(() => { const e = performance.getEntriesByName('first-contentful-paint')[0]; return e ? Math.round(e.startTime) : null })()")
    if plate in ("II",):
        res["molstar_initialised"] = page.evaluate("!!document.querySelector('#molstar canvas')")
    if plate == "I":
        # visible tip labels must not overlap at default zoom (labels hidden below the zoom threshold)
        res["label_overlaps"] = page.evaluate("""(() => {
          const els = [...document.querySelectorAll('#tree text.tip')].filter(t => t.getAttribute('display') !== 'none');
          const boxes = els.map(e => e.getBoundingClientRect()).filter(b => b.width > 0);
          let n = 0;
          for (let i = 0; i < boxes.length; i++) for (let j = i + 1; j < boxes.length; j++) {
            const a = boxes[i], b = boxes[j];
            if (a.left < b.right - 1 && b.left < a.right - 1 && a.top < b.bottom - 1 && b.top < a.bottom - 1) n++;
          }
          return n; })()""")
    shot = OUT / "screens" / f"plate_{plate}_{vp}.png"
    page.screenshot(path=str(shot), full_page=True)
    res["screenshot"] = str(shot.relative_to(OUT.parent))
    if Axe is not None and vp == "desktop":
        r = Axe().run(page)
        viol = [v for v in r.response.get("violations", []) if v.get("impact") in ("serious", "critical")]
        res["axe_serious_critical"] = [{"id": v["id"], "impact": v["impact"], "n": len(v["nodes"]),
                                        "help": v["help"]} for v in viol]
    res["console_errors"] = errors[:10]
    return res


def main() -> int:
    (OUT / "screens").mkdir(parents=True, exist_ok=True)
    results = []
    with sync_playwright() as p:
        browser = p.chromium.launch(channel="chrome", headless=True)
        for vp, (w, h) in VIEWPORTS.items():
            ctx = browser.new_context(viewport={"width": w, "height": h}, device_scale_factor=1)
            for plate, path in PLATES:
                page = ctx.new_page()
                try:
                    results.append(check_page(page, plate, path, vp))
                except Exception as e:  # noqa: BLE001
                    results.append({"plate": plate, "path": path, "viewport": vp, "error": str(e)[:300]})
                page.close()
            ctx.close()
        browser.close()
    fails = []
    for r in results:
        tag = f"{r['plate']}@{r['viewport']}"
        if r.get("error"):
            fails.append(f"{tag}: {r['error']}")
            continue
        if r.get("status") != 200:
            fails.append(f"{tag}: HTTP {r.get('status')}")
        if r.get("horizontal_scroll"):
            fails.append(f"{tag}: horizontal scroll")
        if r.get("console_errors"):
            fails.append(f"{tag}: console errors {r['console_errors'][:2]}")
        if r.get("label_overlaps"):
            fails.append(f"{tag}: {r['label_overlaps']} overlapping tree labels")
        if r.get("molstar_initialised") is False:
            fails.append(f"{tag}: Mol* did not initialise")
        if r.get("axe_serious_critical"):
            fails.append(f"{tag}: axe {[v['id'] for v in r['axe_serious_critical']]}")
        if r.get("fcp_ms") and r["fcp_ms"] > 1500:
            fails.append(f"{tag}: FCP {r['fcp_ms']} ms")
    (OUT / "visual_qc.json").write_text(json.dumps({"results": results, "failures": fails}, indent=1))
    print(json.dumps(fails, indent=1))
    print(f"{len(results)} page checks, {len(fails)} failures")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
