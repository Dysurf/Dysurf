#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
SRC_ROOT="$(cd "$ROOT/../../.." && pwd)"

cd "$SRC_ROOT/src/src_gfortran"
make mpi

cd "$ROOT"
./../../../src/src_gfortran/dysurf PdSe2_CNCS12.txt
mv -f SQE_300K.dat SQE_300K_CNCS12.dat

./../../../src/src_gfortran/dysurf PdSe2_psf2d_hcol_20.txt
mv -f SQE_300K.dat SQE_300K_psf_hcol_20.dat

./../../../src/src_gfortran/dysurf PdSe2_psf2d_hcol_80_120.txt
mv -f SQE_300K.dat SQE_300K_psf_hcol_80_120.dat

source /home/davy/miniconda3/etc/profile.d/conda.sh
conda activate neutronpy
python compare_cncs12_psf_hcol.py
