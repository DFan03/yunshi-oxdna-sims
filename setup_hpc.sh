#!/usr/bin/env bash
# Setup on TU Delft DelftBlue (login node). Run:  bash setup_hpc.sh [--build-oxdna]
# Creates a venv under ENV_ROOT (scratch, not $HOME), installs requirements, smoke-tests.
# --build-oxdna: also builds oxDNA with CUDA + Python bindings at commit 8eed395 into OXDNA_DIR (opt-in; slow, login node).
# No sudo, no credentials. Module names were checked with `module spider` on DelftBlue (see docs/PORTABILITY.md).
set -euo pipefail
cd "$(dirname "$0")"

# ---- configure here (or override from the environment) ----
SCRATCH_BASE=${SCRATCH_BASE:-/scratch/$USER}            # e.g. /scratch/donglifan
ENV_ROOT=${ENV_ROOT:-$SCRATCH_BASE/yunshi-oxdna-env}    # python venv location
OXDNA_DIR=${OXDNA_DIR:-$SCRATCH_BASE/oxDNA}             # oxDNA source + build (only with --build-oxdna)
MODULES_STACK=${MODULES_STACK:-2025}
MODULES=${MODULES:-"gcc/12.4.0 cmake/3.30.5 cuda/12.5 python/3.11.9"}
OXDNA_COMMIT=8eed395357533d0a81aecfb73a0f228aab4f180c
# -----------------------------------------------------------
BUILD=0
for a in "$@"; do [ "$a" = "--build-oxdna" ] && BUILD=1; done

module load "$MODULES_STACK"
# shellcheck disable=SC2086
module load $MODULES
echo "[setup] modules: $(module -t list 2>&1 | tr '\n' ' ')"

if [ ! -x "$ENV_ROOT/bin/python" ]; then
  mkdir -p "$(dirname "$ENV_ROOT")"
  python3 -m venv "$ENV_ROOT"
fi
"$ENV_ROOT/bin/python" -m pip install --upgrade pip
"$ENV_ROOT/bin/python" -m pip install -r requirements.txt

if [ "$BUILD" = 1 ]; then
  echo "[setup] building oxDNA (CUDA) at $OXDNA_COMMIT in $OXDNA_DIR"
  [ -d "$OXDNA_DIR" ] || git clone https://github.com/lorenzo-rovigatti/oxDNA.git "$OXDNA_DIR"
  git -C "$OXDNA_DIR" checkout "$OXDNA_COMMIT"
  mkdir -p "$OXDNA_DIR/build"
  ( cd "$OXDNA_DIR/build" \
    && cmake .. -DCUDA=ON -DCUDA_COMMON_ARCH=ON -DPython=ON -DPython3_EXECUTABLE="$ENV_ROOT/bin/python" \
    && make -j"${SLURM_CPUS_ON_NODE:-4}" && make install )
  echo "[setup] oxDNA binary: $OXDNA_DIR/build/bin/oxDNA  (scripts hard-code other paths: see docs/PORTABILITY.md)"
fi

"$ENV_ROOT/bin/python" tools/smoke_test.py
cat <<MSG
[setup] done. In every shell / job script:
  module load $MODULES_STACK; module load $MODULES
  source $ENV_ROOT/bin/activate
Existing job scripts in jobs/ use 'conda activate oxdna' instead; see docs/PORTABILITY.md.
MSG
