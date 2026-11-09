#!/usr/bin/env python3
"""Smoke test used by setup_local.sh and setup_hpc.sh: imports + py_compile of every script.
Usage: python tools/smoke_test.py   (run from repo root)"""
import importlib, py_compile, sys, pathlib

assert sys.version_info >= (3, 9), "python>=3.9 required"
for m in ["numpy", "scipy", "pandas", "matplotlib", "seaborn"]:
    mod = importlib.import_module(m)
    print(f"  ok {m} {getattr(mod, '__version__', '')}")
for m in ["oxDNA_analysis_tools", "oxpy"]:
    try:
        importlib.import_module(m); print(f"  ok {m}")
    except Exception as e:
        print(f"  WARNING {m} not importable ({type(e).__name__}: {e}); scripts using it will not run "
              "until oxDNA is built with Python bindings (README, 'Installing oxDNA').")
root = pathlib.Path(__file__).resolve().parent.parent
files = [f for d in ("rna", "dna", "tools", "archive") for f in (root / d).rglob("*.py")]
for f in files:
    py_compile.compile(str(f), doraise=True)
for c in root.rglob("__pycache__"):
    import shutil; shutil.rmtree(c, ignore_errors=True)
print(f"  compiled {len(files)} python files OK")
