#!/usr/bin/env bash
# Ensure a single OpenMP runtime: point pip torch's bundled libomp at conda's libomp.
set -euo pipefail
ENV="$(dirname "$0")/../.pixi/envs/default"
T="$ENV/lib/python3.12/site-packages/torch/lib/libomp.dylib"
if [ -f "$T" ] && [ ! -L "$T" ]; then mv "$T" "$T.bundled"; ln -s "$(cd "$ENV/lib" && pwd)/libomp.dylib" "$T"; echo "linked torch libomp to conda libomp"; else echo "already linked"; fi
