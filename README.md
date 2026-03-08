# **Dysurf**, a program for simulating four-dimensional dynamical structure factors

Version: 1.1

Copyright (C) 2023-2025 Yongheng Li <davy_li96@163.com>

Copyright (C) 2019-2023 Changpeng Lin <changpeng.lin@epfl.ch>

Copyright (C) 2019-2025 Bin Wei <binwei@hpu.edu.cn>

Copyright (C) 2019-2025 Jiawang Hong <hongjw@bit.edu.cn>

The Dysurf program is mainly written in Fortran 90 and use some Fortran 2003 extensions.
In this distribution, it contains four subdirectories:

- `/src/src_intel`: Fortran source codes of Dysurf program, this is intel mkl version.

- `/src/src_gfortran`: Fortran source codes of Dysurf program, this is gfortran openblas version. Sometimes it return Errors as below:

```
/usr/bin/ld: /home/davy/software/ds/Dysurf/src/src_gfortran/qpoints.f90:51: undefined reference to `dgetri_'
/usr/bin/ld: /home/davy/software/ds/Dysurf/src/src_gfortran/qpoints.f90:61: undefined reference to `dnrm2_'
/usr/bin/ld: /home/davy/software/ds/Dysurf/src/src_gfortran/qpoints.f90:70: undefined reference to `dnrm2_'
/usr/bin/ld: /home/davy/software/ds/Dysurf/src/src_gfortran/qpoints.f90:104: undefined reference to `dnrm2_'
collect2: error: ld returned 1 exit status
```
Run `gfortran -g -O2 -ffree-line-length-none -fPIE -c dysurf.f90 -o dysurf` seems to resolve this problem. Or you may need to point to your own openblas path.

Choose any version that matches your local system environment.

`/docs`: a user manual for description of running a calculation. In `/docs/Dysurf_online_version_manual` are manual for using Dysurf online version.

`/examples`: examples of running Dysurf program with input files and results

`/tools`: some python scripts to do preprocessing and plot results

Alternatively, try our web interface (http://36.138.185.163:5000/) to bypass compilation problems. The tutorial for using the Dysurf online version can be found in `./docs/Dysurf_online_version_manual`, with corresponding examples located in `examples/CsI/Dysurf_online_version`.

**Note**

Large-scale computations (e.g., high-density q-point sampling or extensive q-space exploration) should not be executed on the web interface.

For computationally intensive tasks, please compile the source code and perform calculations on your local machine.

If you find this code useful, we would appreciate a citation to: XXX (To be published)

Hope you enjoy this program. Please do not hesitate to contact us if any problem rises.


