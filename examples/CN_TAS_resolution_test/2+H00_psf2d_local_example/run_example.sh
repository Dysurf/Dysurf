#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
SRC_ROOT="$(cd "$ROOT/../../.." && pwd)"

cd "$SRC_ROOT/src/src_gfortran"
make mpi

cd "$ROOT"
source /home/davy/miniconda3/etc/profile.d/conda.sh
conda activate neutronpy
python run_psf2d_example.py
