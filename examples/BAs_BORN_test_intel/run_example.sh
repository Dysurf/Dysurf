#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
SRC_ROOT="$(cd "$ROOT/../.." && pwd)"

set +u
source /home/davy/miniconda3/etc/profile.d/conda.sh
conda activate neutronpy
set -u

cd "$SRC_ROOT/src/src_intel"
make mpi

cd "$ROOT"
MPI_LAUNCHER="${MPI_LAUNCHER:-/usr/bin/mpirun}"
"$MPI_LAUNCHER" -np 1 ../../src/src_intel/dysurf BAs_BORN_test.txt
