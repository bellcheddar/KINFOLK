#!/usr/bin/env bash
# Pull-based deploy for KINFOLK (as in PODIUM). Runs on the droplet every minute from
# kinfolk-pull.timer: if the public `live` branch on GitHub has moved, publish its app/ folder.
# The Mac only ever talks HTTPS (git push); nothing inbound over SSH is needed to deploy.
set -euo pipefail

BASE=/opt/kinfolk-deploy
REPO="$BASE/repo"
APP=/srv/kinfolk/app
LOG=/var/log/kinfolk-pull.log
URL=https://github.com/bellcheddar/KINFOLK.git

exec 9>/run/kinfolk-pull.lock
flock -n 9 || exit 0                              # a slow run is never doubled up

mkdir -p "$BASE"
[[ -d "$REPO/.git" ]] || git init -q "$REPO"
cd "$REPO"
git fetch -q --depth 1 "$URL" live               # live is one snapshot, replaced each deploy
new=$(git rev-parse FETCH_HEAD)
cur=$(cat "$BASE/deployed" 2>/dev/null || true)
[[ "$new" == "$cur" ]] && exit 0

git -c advice.detachedHead=false checkout -q --force FETCH_HEAD
src=$(cat source-commit.txt 2>/dev/null || echo "$new")

mkdir -p "$APP"
old_lock=$(sha256sum "$APP/uv.lock" 2>/dev/null | cut -d' ' -f1 || true)
# Standard modes whatever the committer's machine had (the GOBSMACKED lesson); keep the venv.
rsync -a --delete --chmod=D755,F644 --exclude '.venv/' app/ "$APP/"
printf '%s\n' "$src" > "$APP/static/version.txt"
new_lock=$(sha256sum "$APP/uv.lock" | cut -d' ' -f1)

if [[ ! -x "$APP/.venv/bin/gunicorn" || "$old_lock" != "$new_lock" ]]; then
  # System Python, so the venv never points into /root (unreadable by the service user).
  (cd "$APP" && UV_PYTHON_PREFERENCE=only-system /usr/local/bin/uv sync --no-dev --frozen -q)
fi

chown -R kinfolk:www-data /srv/kinfolk
systemctl restart kinfolk.service
printf '%s\n' "$new" > "$BASE/deployed"
echo "$(date -Is) deployed live=$new source=$src" >> "$LOG"
