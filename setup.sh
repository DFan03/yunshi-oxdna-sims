#!/usr/bin/env bash
# One-click setup: creates env "yunshi-oxdna" (conda/mamba) or ./.venv (venv+pip), installs deps,
# runs a smoke test (imports + py_compile on every script).
# oxDNA itself (the simulator) is an EXTERNAL prerequisite: see README "Installing oxDNA".
set -euo pipefail
cd "$(dirname "$0")"
ENV_NAME=yunshi-oxdna

if command -v mamba >/dev/null 2>&1 || command -v conda >/dev/null 2>&1; then
  TOOL=$(command -v mamba || command -v conda)
  echo "[setup] using $TOOL"
  "$TOOL" env create -n "$ENV_NAME" -f environment.yml -y || "$TOOL" env update -n "$ENV_NAME" -f environment.yml
  PY=$("$TOOL" run -n "$ENV_NAME" python -c 'import sys;print(sys.executable)')
  ACTIVATE="conda activate $ENV_NAME"
else
  echo "[setup] conda/mamba not found, using venv + pip"
  PYBIN=${PYTHON:-python3}
  "$PYBIN" -m venv .venv
  PY="$PWD/.venv/bin/python"
  "$PY" -m pip install --upgrade pip
  "$PY" -m pip install -r requirements.txt
  ACTIVATE="source .venv/bin/activate"
fi

echo "[setup] smoke test: imports"
"$PY" - <<'PYEOF'
import importlib, sys
assert sys.version_info >= (3, 9), "python>=3.9 required"
core = ["numpy", "scipy", "pandas", "matplotlib", "seaborn"]
for m in core:
    mod = importlib.import_module(m); print(f"  ok {m} {getattr(mod,'__version__','')}")
# oxDNA tools: need the oxpy bindings (shipped with an oxDNA build or the 'oxpy' wheel)
for m in ["oxDNA_analysis_tools", "oxpy"]:
    try:
        importlib.import_module(m); print(f"  ok {m}")
    except Exception as e:
        print(f"  WARNING {m} not importable ({type(e).__name__}: {e}); analysis scripts that use it will not run "
              "until oxDNA is built with Python bindings (see README).")
PYEOF

echo "[setup] smoke test: py_compile on all scripts"
find rna dna tools archive -name '*.py' -print0 | xargs -0 "$PY" -m py_compile
find . -name '__pycache__' -type d -prune -exec rm -rf {} + 2>/dev/null || true
n=$(find rna dna tools archive -name '*.py' | wc -l | tr -d ' ')
echo "[setup] compiled $n python files OK"
echo "[setup] done. Activate with:  $ACTIVATE"
