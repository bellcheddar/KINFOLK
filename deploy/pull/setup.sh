#!/usr/bin/env bash
# One-time droplet setup for KINFOLK's HTTPS pull deploys. Idempotent: rerun it freely.
# Run as root from the DigitalOcean web console with one typed line:
#   curl -fsSL https://raw.githubusercontent.com/bellcheddar/KINFOLK/main/deploy/pull/setup.sh | bash
# The repo is public, so the droplet pulls over plain HTTPS: no deploy key, no inbound SSH.
set -euo pipefail
RAW=https://raw.githubusercontent.com/bellcheddar/KINFOLK/main/deploy
SITE=kinfolk.mdeller.com

echo "==> Prerequisites"
for b in git rsync curl nginx systemctl; do command -v "$b" >/dev/null || { echo "missing: $b"; exit 1; }; done
id kinfolk >/dev/null 2>&1 || useradd --system --no-create-home --shell /usr/sbin/nologin kinfolk
if [[ ! -x /usr/local/bin/uv ]]; then
  curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin INSTALLER_NO_MODIFY_PATH=1 sh
fi
python3 -c 'import sys; assert sys.version_info >= (3, 12), sys.version' || { echo "KINFOLK needs Python >= 3.12"; exit 1; }

echo "==> Service, nginx site, pull script and timer"
curl -fsSL "$RAW/kinfolk.service"            -o /etc/systemd/system/kinfolk.service
curl -fsSL "$RAW/pull/kinfolk-pull.service"  -o /etc/systemd/system/kinfolk-pull.service
curl -fsSL "$RAW/pull/kinfolk-pull.timer"    -o /etc/systemd/system/kinfolk-pull.timer
curl -fsSL "$RAW/pull/kinfolk-pull.sh"       -o /usr/local/bin/kinfolk-pull
chmod 755 /usr/local/bin/kinfolk-pull
if [[ ! -f /etc/nginx/sites-available/$SITE ]] || ! grep -q "ssl_certificate" "/etc/nginx/sites-available/$SITE"; then
  curl -fsSL "$RAW/nginx-kinfolk.conf" -o "/etc/nginx/sites-available/$SITE"   # certbot adds TLS below
fi
ln -sf "/etc/nginx/sites-available/$SITE" "/etc/nginx/sites-enabled/$SITE"
# Certbot owns the TLS lines in that file, so it is never re-downloaded once it has them.
# Move any older TCP upstream onto the socket in place (the droplet runs several apps, and
# 127.0.0.1:8010 belongs to one of them).
sed -i 's|proxy_pass http://127\.0\.0\.1:8010;|proxy_pass http://unix:/run/kinfolk/kinfolk.sock;|' \
  "/etc/nginx/sites-available/$SITE"
systemctl daemon-reload
systemctl enable kinfolk.service

echo "==> First publish of the live branch (installs the Python environment; a few minutes)"
rm -f /opt/kinfolk-deploy/deployed
/usr/local/bin/kinfolk-pull

echo "==> nginx and TLS"
nginx -t
systemctl reload nginx
if [[ ! -d /etc/letsencrypt/live/$SITE ]]; then
  certbot --nginx -d "$SITE" --redirect -m marc@marcdeller.com --agree-tos -n
fi
systemctl enable --now kinfolk-pull.timer

echo "==> Checks"
for _ in $(seq 1 30); do
  systemctl is-active --quiet kinfolk.service && break
  sleep 2
done
if ! systemctl is-active --quiet kinfolk.service; then
  echo "kinfolk.service did not start:"
  journalctl -u kinfolk.service -n 30 --no-pager
  exit 1
fi
for p in / /k/ABL1 /almanac /about /static/version.txt; do
  printf '  %-22s %s\n' "$p" "$(curl -s -o /dev/null -w '%{http_code}' "https://$SITE$p")"
done
echo "KINFOLK_READY live=$(cat /opt/kinfolk-deploy/deployed) source=$(cat /srv/kinfolk/app/static/version.txt)"
