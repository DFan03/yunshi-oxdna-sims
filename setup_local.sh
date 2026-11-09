#!/usr/bin/env bash
# Setup for a laptop/workstation (macOS/Linux).
# Creates conda/mamba env "yunshi-oxdna" if conda/mamba exist, else ./.venv (venv+pip);
# installs requirements; runs a smoke test.
# Option: --build-oxdna   also builds oxDNA (CPU, with Python bindings) at commit 8eed395 into ./oxDNA
#                         (needs git, cmake, a C++ compiler; untested by the authors on every platform).
set -euo pipefail
cd "$(dirname "$0")"
ENV_NAME=yunshi-oxdna
OXDNA_COMMIT=8eed395357533d0a81aecfb73a0f228aab4f180c
BUILD=0
for a in "$@"; do [ "$a" = "--build-oxdna" ] && BUILD=1; done

if command -v mamba >/dev/null 2>&1 || command -v conda >/dev/null 2>&1; then
  TOOL=$(command -v mamba || command -v conda)
  echo "[setup] using $TOOL"
  "$TOOL" env create -n "$ENV_NAME" -f environment.yml -y || "$TOOL" env update -n "$ENV_NAME" -f environment.yml
  RUN=("$TOOL" run -n "$ENV_NAME")
  ACTIVATE="conda activate $ENV_NAME"
else
  echo "[setup] conda/mamba not found, using venv + pip"
  python3 -m venv .venv
  .venv/bin/python -m pip install --upgrade pip
  .venv/bin/python -m pip install -r requirements.txt
  RUN=("$PWD/.venv/bin/python")
  ACTIVATE="source .venv/bin/activate"
fi

# run python through the chosen environment
if [ "${#RUN[@]}" = 1 ]; then PYRUN=("${RUN[@]}"); else PYRUN=("${RUN[@]}" python); fi

if [ "$BUILD" = 1 ]; then
  echo "[setup] building oxDNA (CPU) at $OXDNA_COMMIT"
  [ -d oxDNA ] || git clone https://github.com/lorenzo-rovigatti/oxDNA.git oxDNA
  git -C oxDNA checkout "$OXDNA_COMMIT"
  mkdir -p oxDNA/build
  PYBIN=$("${PYRUN[@]}" -c 'import sys;print(sys.executable)')
  ( cd oxDNA/build && cmake .. -DPython=ON -DPython3_EXECUTABLE="$PYBIN" && make -j"$(getconf _NPROCESSORS_ONLN)" && make install )
  echo "[setup] oxDNA binary: $PWD/oxDNA/build/bin/oxDNA"
fi

echo "[setup] smoke test"
"${PYRUN[@]}" tools/smoke_test.py
echo "[setup] done. Activate with:  $ACTIVATE"
