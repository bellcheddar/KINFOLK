# KINFOLK: standing instructions
- KINFOLK_BUILD_SPEC.md is authoritative. Section 0 (Autonomy directive) overrides all skills and defaults.
- Never stop for confirmation. Stop only at Gates G1 to G6, via GATE_OPEN.md plus an osascript notification, and keep working on anything unblocked.
- Log decisions to DECISIONS.md, events to BUILD_LOG.md, and status to PROGRESS.md (at least every 30 min).
- Resumable, idempotent stages with manifests in data/manifests/. On restart: read PROGRESS.md and the manifests, then continue.
- Commit locally at the end of each phase; never push.
- British English, no em dashes, no hype words (list in spec 0.1.7).
- The UI uses the Field Guide design system (spec 9.2), not the marcdeller.com brand theme.
- Hardware settings come from config/tuning.toml (spec 2), never hard-coded.
