#!/usr/bin/env bash
# Emergency repair: restore ownership and standard modes on the droplet, restart, print page codes.
#   bash deploy/repair_perms.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
if [[ -f .env ]]; then set -a; source .env; set +a; fi
DROPLET_SSH="${DROPLET_SSH:-root@mdeller.com}"
# One multiplexed connection for every step: the droplet rate-limits and bans IPs that open
# several SSH connections in quick succession (a 5-connection deploy got banned, 2026-09-28).
CTL="$HOME/.ssh/cm-kinfolk-%C"
SSH=(ssh -o BatchMode=yes -o ConnectTimeout=15 -o ControlMaster=auto -o ControlPath="$CTL" -o ControlPersist=120)
[[ -n "${SSH_KEY:-}" ]] && SSH+=(-i "${SSH_KEY/#\~/$HOME}")
"${SSH[@]}" "$DROPLET_SSH" bash -s <<'REMOTE_EOF'
set -euo pipefail
D=/srv/kinfolk
find "$D" -path "$D/app/.venv" -prune -o -exec chown kinfolk:www-data {} +
find "$D" -path "$D/app/.venv" -prune -o -name .env -prune -o -type d -exec chmod 755 {} + -o -type f -exec chmod 644 {} +
systemctl restart kinfolk.service
sleep 2
systemctl is-active kinfolk.service
for p in / /k/ABL1 /almanac /about; do printf '%-10s %s\n' "$p" "$(curl -s -o /dev/null -w '%{http_code}' http://127.0.0.1:8010$p)"; done
REMOTE_EOF
