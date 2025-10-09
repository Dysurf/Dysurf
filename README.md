### **Dysurf**, a program for simulating four-dimensional dynamical structure factors

Version: 1.1

Copyright (C) 2023-2025 Yongheng Li <davy_li96@163.com>

Copyright (C) 2019-2023 Changpeng Lin <changpeng.lin@epfl.ch>

Copyright (C) 2019-2023 Jiawang Hong <hongjw@bit.edu.cn>

The Dysurf program is mainly written in Fortran 90 and use some Fortran 2003 extensions.
In this distribution, it contains four subdirectories:

- `/src_intel`: Fortran source codes of Dysurf program, this is intel mkl version.

- `/src_gfortran`: Fortran source codes of Dysurf program, this is gfortran openblas version. 

Choose any version that matches your local system environment.

`/docs`: a user manual for description of running a calculation

`/examples`: examples of running Dysurf program with input files and results

`/tools`: some python scripts to do preprocessing and plot results

Alternatively, try our web interface (http://36.138.185.163:5000/) to bypass compilation problems.

Note:

Large-scale computations (e.g., high-density q-point sampling or extensive q-space exploration) should not be executed on the web interface.

For computationally intensive tasks, please compile the source code and perform calculations on your local machine.

If you find this code useful, we would appreciate a citation to: XXX (To be published)

Hope you enjoy this program. Please do not hesitate to contact us if any problem rises.


