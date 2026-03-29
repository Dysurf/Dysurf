# **Dysurf**, a program for simulating four-dimensional dynamical structure factors

Version: 1.2

Copyright (C) 2023-2025 Yongheng Li <davy_li96@163.com>

Copyright (C) 2019-2023 Changpeng Lin <changpeng.lin@epfl.ch>

Copyright (C) 2019-2025 Bin Wei <binwei@hpu.edu.cn>

Copyright (C) 2019-2025 Jiawang Hong <hongjw@bit.edu.cn>

The Dysurf program is mainly written in Fortran 90 and use some Fortran 2003 extensions.
In this distribution, it contains four subdirectories:

- `/src/src_intel`: Fortran source codes of Dysurf program, this is intel mkl version.

- `/src/src_gfortran`: Fortran source codes of Dysurf program, this is gfortran openblas version. Sometimes it return Errors as below:

If `OpenMPI` is available on your machine, you can also compile the MPI-enabled gfortran version in `/src/src_gfortran` and run Dysurf with `mpirun`.

Choose any version that matches your local system environment.

## Compilation

The currently validated build paths are:

### 1. gfortran + OpenBLAS

```bash
cd your_own_path/Dysurf/src/src_gfortran
make
```

This builds the serial executables:

- `dysurf`
- `validate_psf2d_parametrization`

### 2. gfortran + OpenMPI

If `mpif90` and `mpirun` are available:

```bash
cd your_own_path/Dysurf/src/src_gfortran
make mpi
```

Then examples can be run with commands such as:

```bash
mpirun -np 1 ./dysurf input.txt
```

### 3. Intel ifx + Intel MPI in conda

One workable installation route is:

```bash
conda install -c https://software.repos.intel.com/python/conda/ -c conda-forge impi-devel
conda install -c https://software.repos.intel.com/python/conda/ -c conda-forge ifx_linux-64
```

If you use a dedicated environment, for example `dysurf`, activate it first:

```bash
source /home/davy/miniconda3/etc/profile.d/conda.sh
conda activate dysurf
```

Then build the Intel version:

```bash
cd your_own_path/Dysurf/src/src_intel
make mpi
```

This now builds:

- `dysurf`
- `validate_psf2d_parametrization`

Notes:

- The current `src/src_intel` `Makefile` is configured to use `mpiifx` when `make mpi` is called.
- If MKL is not present in the active environment, the Intel build currently falls back to system `lapack/openblas` for linking.
- Running Intel-MPI executables should use the matching Intel MPI launcher from the same environment, rather than mixing them with a system OpenMPI `mpirun`.

`/docs`: a user manual for description of running a calculation. In `/docs/Dysurf_online_version_manual` are manual for using Dysurf online version.

`/examples`: examples of running Dysurf program with input files and results

`/tools`: some python scripts to do preprocessing and plot results

Alternatively, try our web interface (http://36.138.185.163:5000/) to bypass compilation problems. The tutorial for using the Dysurf online version can be found in `./docs/Dysurf_online_version_manual`, with corresponding examples located in `examples/CsI/Dysurf_online_version`.

**Note**

Large-scale computations (e.g., high-density q-point sampling or extensive q-space exploration) should not be executed on the web interface.

For computationally intensive tasks, please compile the source code and perform calculations on your local machine.

If you find this code useful, we would appreciate a citation to: XXX (To be published)

Hope you enjoy this program. Please do not hesitate to contact us if any problem rises.
