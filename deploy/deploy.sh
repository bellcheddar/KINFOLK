#!/usr/bin/env bash
# Deploy KINFOLK over HTTPS (as in PODIUM). Run from anywhere: bash deploy/deploy.sh
#
#   1. tests pass locally
#   2. the committed code goes to GitHub `main`; app/ plus the data bundle (kept out of main)
#      goes to `live` as a single snapshot commit, force-replaced each deploy so the repo
#      does not grow by 150 MB per release
#   3. the droplet's kinfolk-pull timer (every minute) publishes it: nothing inbound over SSH
#   4. this script waits until https://kinfolk.mdeller.com/static/version.txt reports the
#      commit, then checks every page (deploy/check_live.sh)
#
# One-time droplet setup: deploy/pull/README.md. Emergency fallback: deploy/deploy-rsync.sh.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
SERVER_NAME="${SERVER_NAME:-kinfolk.mdeller.com}"
WAIT=1
[[ "${1:-}" == "--no-wait" ]] && WAIT=0   # first deploy: push, then run deploy/pull/setup.sh on the droplet

echo "==> Tests"
pixi run pytest pipeline/tests -q
(cd app && uv run pytest tests -q)

if [[ -n "$(git status --porcelain -- app pipeline deploy)" ]]; then
  echo "Uncommitted changes in app/, pipeline/ or deploy/: commit them first (a deploy publishes a commit)."
  exit 1
fi
[[ -f app/data/kinfolk.sqlite ]] || { echo "No data bundle: run kinfolk run --phase 8 first."; exit 1; }
sha=$(git rev-parse HEAD)

echo "==> Pushing ${sha:0:12} to origin/main (HTTPS)"
git push -q origin HEAD:refs/heads/main

echo "==> Building the live snapshot (app/ + data bundle)"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
export GIT_INDEX_FILE="$tmp/index"
git add -f -- app ':!app/.venv' ':!app/tests' ':!app/__pycache__' ':!(glob)app/**/__pycache__/**' ':!(glob)app/**/.DS_Store'
blob=$(printf '%s\n' "$sha" | git hash-object -w --stdin)
git update-index --add --cacheinfo "100644,$blob,source-commit.txt"
tree=$(git write-tree)
unset GIT_INDEX_FILE
live=$(git commit-tree "$tree" -m "live: app and data bundle for ${sha:0:12}")
echo "    $(git ls-tree -r --name-only "$live" | wc -l | tr -d ' ') files, snapshot ${live:0:12}"

echo "==> Pushing the snapshot to origin/live"
git push -q --force origin "$live:refs/heads/live"

if [[ $WAIT -eq 0 ]]; then
  echo "==> Pushed. Now run the one-time setup on the droplet (deploy/pull/README.md)."
  exit 0
fi

echo "==> Waiting for the droplet to publish it (it checks every minute; first run installs Python)"
live_sha=""
for i in $(seq 1 60); do
  live_sha=$(curl -s --max-time 10 "https://${SERVER_NAME}/static/version.txt?t=$(date +%s)" | tr -d '[:space:]' || true)
  if [[ "$live_sha" == "$sha" ]]; then echo "    live after ~$((i * 10)) s"; break; fi
  sleep 10
done
[[ "$live_sha" == "$sha" ]] || { echo "    not live after 10 min: see /var/log/kinfolk-pull.log on the droplet"; exit 1; }

echo "==> Checking every page"
bash deploy/check_live.sh "https://${SERVER_NAME}"
echo "==> Live: https://${SERVER_NAME} at ${sha:0:12}"
