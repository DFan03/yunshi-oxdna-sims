#!/usr/bin/env bash
# Thin wrapper: DelftBlue (sbatch present or hostname contains 'delftblue') -> setup_hpc.sh, else setup_local.sh.
cd "$(dirname "$0")"
if command -v sbatch >/dev/null 2>&1 || [[ "$(hostname -f 2>/dev/null || hostname)" == *delftblue* ]]; then
  exec bash setup_hpc.sh "$@"
fi
exec bash setup_local.sh "$@"
