# Gate G3 open: deploy approval

## What is needed from Marc
Approval to deploy KINFOLK to the droplet (mdeller.com, 45.55.102.228) as `kinfolk.mdeller.com`: app and 0.16 GB data bundle to `/srv/kinfolk`, a `kinfolk` systemd service (gunicorn on 127.0.0.1:8010), an nginx site, then DNS and TLS.

## How to unblock
1. Approve the deploy scripts for Claude Code via `/permissions` (allow rules):
   - `Bash(bash deploy/deploy.sh)`
   - `Bash(bash deploy/repair_perms.sh)`
   - `Bash(bash deploy/check_live.sh)`
2. Say "deploy KINFOLK". Claude runs `bash deploy/deploy.sh`, which first checks port 22, the known_hosts entry and one SSH login (stopping without retrying if any fail), syncs with standard file modes, fixes ownership before installing, starts the service and checks every page on the droplet.
3. Add a DNS A record: `kinfolk.mdeller.com -> 45.55.102.228`.
4. Issue the certificate (the deploy script prints the exact certbot command), then Claude runs `bash deploy/check_live.sh`.

SSH uses the `~/.ssh/config` entry for mdeller.com (key `~/.ssh/faffabout_deploy`); nothing else is needed in `.env` (see `.env.example`).

## Completed
- Phases 0 to 10: 1,325 of 1,326 fold batches (one failed: MASTL, SVD did not converge); annotation, pockets, selectivity, induced fit, findings, bundle.
- QC: 20 pipeline and 16 app tests pass; visual QC 16/16 pages at desktop and mobile width; axe clean; science checks all pass (`qc/science.md`).
- Headline results in `FINDINGS.md`.

## Before approving, note
- The droplet also runs GOBSMACKED (up to 1.4 GB RAM of 3.8 GB); KINFOLK's two gunicorn workers need roughly 200 to 300 MB.
- Deployment writes system files as root (`/etc/systemd/system/kinfolk.service`, `/etc/nginx/sites-available/kinfolk.mdeller.com`) and reloads nginx, which serves the other sites too; `nginx -t` runs before the reload.

## After the gate
Only G4 (publishing to GitHub) remains, if wanted.
