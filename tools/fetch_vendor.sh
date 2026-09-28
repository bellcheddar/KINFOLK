#!/usr/bin/env bash
# Fetch third-party tools used by the pipeline into tools/vendor/ (not committed).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p vendor
[ -d vendor/p2rank_2.5.1 ] || curl -sL https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz | tar -xz -C vendor
[ -d vendor/Kincore-standalone ] || git clone -q --depth 1 https://github.com/vivekmodi/Kincore-standalone vendor/Kincore-standalone
echo "vendor tools ready"
