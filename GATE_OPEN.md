# Gate G3 open: deploy approval (HTTPS pull deploys, as in PODIUM)

## What is needed from Marc
1. Create the public repo (G4, approved as public): `gh repo create bellcheddar/KINFOLK --public`.
2. After Claude's first push, type one line in the DigitalOcean web console (as root):
   `curl -fsSL https://raw.githubusercontent.com/bellcheddar/KINFOLK/main/deploy/pull/setup.sh | bash`

## How it works
- `bash deploy/deploy.sh` (Mac): tests, pushes code to `main` and a snapshot of `app/` plus the
  data bundle to `live` over HTTPS, waits for `/static/version.txt`, checks every page.
- The droplet's `kinfolk-pull.timer` fetches `live` every minute and publishes it. No inbound SSH.
- First time only: `bash deploy/deploy.sh --no-wait`, then the console line above (it installs the
  service, nginx site, pull timer and TLS certificate, and ends `KINFOLK_READY`).
- Details and operations: `deploy/pull/README.md`. SSH fallback: `deploy/deploy-rsync.sh`.

## Completed
- Phases 0 to 10: 1,325 of 1,326 fold batches; annotation, pockets, selectivity, induced fit,
  findings, bundle (164 MB snapshot, 3,254 files).
- QC: 20 pipeline and 16 app tests; visual QC 16/16; axe clean; science checks pass.
- Pre-publication scan: no secrets or personal network details in the history.

## After the gate
Nothing: G3 and G4 close together once https://kinfolk.mdeller.com is live.
