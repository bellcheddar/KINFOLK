#!/usr/bin/env bash
# Deploy KINFOLK to the droplet. NOT run by the build: this is gate G3.
#
#   bash deploy/deploy.sh
#
# Reads DROPLET_SSH (default root@mdeller.com) from .env. The SSH key comes from
# ~/.ssh/config (Host mdeller.com 45.55.102.228) or SSH_KEY in .env.
# Follows the vibe-coding Deploy rules: standard file modes on rsync, ownership
# fixed before anything that can fail, pre-flight SSH checks, page-level checks.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if [[ -f .env ]]; then set -a; source .env; set +a; fi
DROPLET_SSH="${DROPLET_SSH:-root@mdeller.com}"
HOST="${DROPLET_SSH#*@}"
REMOTE=/srv/kinfolk
# One multiplexed connection for every step: the droplet rate-limits and bans IPs that open
# several SSH connections in quick succession (a 5-connection deploy got banned, 2026-09-28).
CTL="$HOME/.ssh/cm-kinfolk-%C"
SSH=(ssh -o BatchMode=yes -o ConnectTimeout=15 -o ControlMaster=auto -o ControlPath="$CTL" -o ControlPersist=120)
[[ -n "${SSH_KEY:-}" ]] && SSH+=(-i "${SSH_KEY/#\~/$HOME}")
RSYNC_SSH="${SSH[*]}"

echo "==> Pre-flight"
ssh-keygen -F "$HOST" >/dev/null || { echo "No known_hosts entry for $HOST. Connect once interactively to verify the key. Stopping."; exit 1; }
"${SSH[@]}" "$DROPLET_SSH" true || { echo "SSH login to $DROPLET_SSH failed (unreachable, banned or refused). Not retrying (ban risk). Stopping."; exit 1; }
trap '"${SSH[@]}" -O exit "$DROPLET_SSH" >/dev/null 2>&1 || true' EXIT

echo "==> Syncing app code and data bundle to $DROPLET_SSH:$REMOTE/app"
"${SSH[@]}" "$DROPLET_SSH" "mkdir -p $REMOTE/app"
rsync -az --delete -e "$RSYNC_SSH" --chmod=Du=rwx,Dgo=rx,Fu=rw,Fgo=r \
  --exclude '.venv/' --exclude '__pycache__/' --exclude 'tests/' --exclude '.env' \
  "$ROOT/app/" "$DROPLET_SSH:$REMOTE/app/"
rsync -az -e "$RSYNC_SSH" --chmod=Fu=rw,Fgo=r "$ROOT/deploy/kinfolk.service" "$ROOT/deploy/nginx-kinfolk.conf" \
  "$DROPLET_SSH:/tmp/"

echo "==> Remote setup: ownership first, then venv, systemd and nginx"
"${SSH[@]}" "$DROPLET_SSH" bash -s <<'REMOTE_EOF'
set -euo pipefail
id kinfolk >/dev/null 2>&1 || useradd --system --no-create-home --shell /usr/sbin/nologin kinfolk
find /srv/kinfolk -path /srv/kinfolk/app/.venv -prune -o -exec chown kinfolk:www-data {} +
cd /srv/kinfolk/app
command -v uv >/dev/null || [ -x "$HOME/.local/bin/uv" ] || curl -LsSf https://astral.sh/uv/install.sh | sh
UV="$(command -v uv || echo "$HOME/.local/bin/uv")"
"$UV" sync --no-dev --frozen 2>/dev/null || "$UV" sync --no-dev
chown -R kinfolk:www-data /srv/kinfolk/app/.venv
mv /tmp/kinfolk.service /etc/systemd/system/kinfolk.service
mv /tmp/nginx-kinfolk.conf /etc/nginx/sites-available/kinfolk.mdeller.com
ln -sf /etc/nginx/sites-available/kinfolk.mdeller.com /etc/nginx/sites-enabled/kinfolk.mdeller.com
systemctl daemon-reload && systemctl enable --now kinfolk.service && systemctl restart kinfolk.service
nginx -t && systemctl reload nginx
sleep 2
systemctl is-active kinfolk.service
echo "--- pages via gunicorn on the droplet:"
fail=0
for p in / /k/ABL1 /resemblance /rooms /black-sheep /dark /almanac /about /api/kinases /api/summary/induced_fit; do
  code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 20 "http://127.0.0.1:8010$p")
  printf '%-26s %s\n' "$p" "$code"; [ "$code" = 200 ] || fail=1
done
exit $fail
REMOTE_EOF

echo "==> Done on the droplet. Remaining one-off steps:"
echo "   1. DNS: A record kinfolk.mdeller.com -> $(dig +short "$HOST" | head -1)"
echo "   2. TLS: ${SSH[*]} $DROPLET_SSH 'certbot --nginx -d kinfolk.mdeller.com --redirect -m marc@marcdeller.com --agree-tos -n'"
echo "   3. Check: bash deploy/check_live.sh"
