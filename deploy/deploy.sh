#!/usr/bin/env bash
# Deploy KINFOLK to the droplet. NOT run by the build: this is gate G3.
# Usage: DROPLET=user@host bash deploy/deploy.sh
set -euo pipefail
: "${DROPLET:?set DROPLET=user@host}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE=/srv/kinfolk

echo "1. Sync app code and data bundle"
rsync -az --delete --chmod=Du=rwx,Dgo=rx,Fu=rw,Fgo=r --exclude '.venv' --exclude '__pycache__' --exclude 'tests' "$ROOT/app/" "$DROPLET:$REMOTE/app/"
rsync -az --chmod=Fu=rw,Fgo=r "$ROOT/deploy/kinfolk.service" "$ROOT/deploy/nginx-kinfolk.conf" "$DROPLET:/tmp/"

echo "2. Remote setup (user, venv, systemd, nginx)"
ssh "$DROPLET" bash -s <<'REMOTE_EOF'
set -euo pipefail
id kinfolk >/dev/null 2>&1 || sudo useradd --system --no-create-home --shell /usr/sbin/nologin kinfolk
cd /srv/kinfolk/app
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh
~/.local/bin/uv sync --no-dev || uv sync --no-dev
sudo chown -R kinfolk:www-data /srv/kinfolk
sudo mv /tmp/kinfolk.service /etc/systemd/system/kinfolk.service
sudo mv /tmp/nginx-kinfolk.conf /etc/nginx/sites-available/kinfolk.mdeller.com
sudo ln -sf /etc/nginx/sites-available/kinfolk.mdeller.com /etc/nginx/sites-enabled/kinfolk.mdeller.com
sudo systemctl daemon-reload && sudo systemctl enable --now kinfolk.service && sudo systemctl restart kinfolk.service
sudo nginx -t && sudo systemctl reload nginx
REMOTE_EOF

echo "3. DNS: add an A record kinfolk.mdeller.com -> droplet IP (at the DNS provider), then:"
echo "   ssh $DROPLET 'sudo certbot --nginx -d kinfolk.mdeller.com --redirect -m marc@marcdeller.com --agree-tos -n'"
echo "4. Check: curl -sI https://kinfolk.mdeller.com | head -1"
