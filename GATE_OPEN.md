# Gate G3 open

## What is needed from Marc
Approval to deploy kinfolk.mdeller.com (rsync bundle to the droplet, systemd, nginx, certbot).

## How to unblock
Review deploy/deploy.sh, then run: bash deploy/deploy.sh (it needs SSH access, which the build may not use).

## Completed
Phases 0 to 10 complete; QC in qc/; findings in FINDINGS.md.

## After the gate
Nothing else is pending locally; after deployment, check https://kinfolk.mdeller.com.
