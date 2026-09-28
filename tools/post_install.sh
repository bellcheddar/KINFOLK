#!/usr/bin/env bash
# Run after every `pixi install` / `pixi add`:
# 1. a single OpenMP runtime: point pip torch's bundled libomp at conda's libomp;
# 2. clear the macOS hidden flag on .pth files (Python >= 3.12.14 skips hidden .pth files,
#    which silently breaks the editable kinfolk install).
set -euo pipefail
ENV="$(cd "$(dirname "$0")/.." && pwd)/.pixi/envs/default"
SP="$ENV/lib/python3.12/site-packages"
T="$SP/torch/lib/libomp.dylib"
if [ -f "$T" ] && [ ! -L "$T" ]; then mv "$T" "$T.bundled"; ln -s "$ENV/lib/libomp.dylib" "$T"; echo "linked torch libomp to conda libomp"; fi
chflags nohidden "$SP"/*.pth 2>/dev/null || true
