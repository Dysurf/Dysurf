.SUFFIXES: .f90

#-----------------------------------------------
# Makefile for the program 'Dysurf'.
# Intel compiler version, synced with src_gfortran features.
#-----------------------------------------------

USE_MPI ?= 0
I_MPI_ROOT ?= $(CONDA_PREFIX)
ifeq ($(USE_MPI),1)
FC = mpiifx
else
FC = ifx
endif
FCL = ${FC}

MPI_INC = -I$(I_MPI_ROOT)/include/mpi -I$(I_MPI_ROOT)/include
MPI_LIBS = -L$(I_MPI_ROOT)/lib/release -L$(I_MPI_ROOT)/lib -lmpifort -lmpi -ldl -lrt -lpthread

FFLAGS = -traceback -debug -O2 -static-intel -fpp $(MPI_INC)

MKL_RT_CANDIDATE := $(firstword $(wildcard $(CONDA_PREFIX)/lib/libmkl_rt.so) $(wildcard /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_rt.so))
ifeq ($(MKL_RT_CANDIDATE),)
BLAS_LDFLAGS = -L/usr/lib/x86_64-linux-gnu/ -llapack -lopenblas -lm -lpthread
else
BLAS_LDFLAGS = -L$(dir $(MKL_RT_CANDIDATE)) -qmkl
endif

LDFLAGS = $(BLAS_LDFLAGS) $(MPI_LIBS)

PROGRAM_SRCS = dysurf.f90 validate_psf2d_parametrization.f90
FCSRC = $(filter-out $(PROGRAM_SRCS),$(wildcard *.f90))

MODS = ${FCSRC:.f90=.mod}

OBJS = ${FCSRC:.f90=.o}

default: dysurf validate_psf2d_parametrization

all: dysurf validate_psf2d_parametrization

mpi:
	$(MAKE) USE_MPI=1 dysurf validate_psf2d_parametrization

dysurf: ${OBJS} dysurf.o
	${FCL} -o $@ $^ ${LDFLAGS}

validate_psf2d_parametrization: ${OBJS} validate_psf2d_parametrization.o
	${FCL} -o $@ $^ ${LDFLAGS}

phonon_spectra.o phonon_spectra.mod: variables.mod constants.mod readin.mod func.mod
qpoints.o qpoints.mod: variables.mod
readin.o readin.mod: variables.mod constants.mod func.mod psf_2d_state.mod
rmsd.o rmsd.mod: phonon_spectra.mod func.mod constants.mod variables.mod
material_sqw.o material_sqw.mod: phonon_spectra.mod func.mod constants.mod variables.mod
mpi_helper.o mpi_helper.mod:
psf_2d_state.o psf_2d_state.mod:
psf_2d_projection.o psf_2d_projection.mod: tas_resolution_cn.mod
psf_2d_param_estimation.o psf_2d_param_estimation.mod:
psf_2d_kernel.o psf_2d_kernel.mod:
psf_2d_convolution.o psf_2d_convolution.mod: psf_2d_kernel.mod
psf_2d_resolution_driver.o psf_2d_resolution_driver.mod: tas_resolution_cn.mod psf_2d_projection.mod psf_2d_param_estimation.mod psf_2d_convolution.mod
tas_resolution_sampling.o tas_resolution_sampling.mod: material_sqw.mod variables.mod
tas_local_linearization.o tas_local_linearization.mod: tas_resolution_cn.mod
tas_numeric4d_convolution.o tas_numeric4d_convolution.mod:
tas_observation_operator.o tas_observation_operator.mod: func.mod tas_resolution_cn.mod
sqecalc.o sqecal.mod: func.mod variables.mod constants.mod tas_resolution_cn.mod tas_local_linearization.mod tas_observation_operator.mod mpi_helper.mod psf_2d_state.mod psf_2d_resolution_driver.mod
dysurf.o: constants.mod func.mod variables.mod phonon_spectra.mod qpoints.mod readin.mod rmsd.mod sqecalc.mod mpi_helper.mod psf_2d_state.mod
validate_psf2d_parametrization.o: variables.mod tas_resolution_cn.mod psf_2d_projection.mod psf_2d_param_estimation.mod psf_2d_kernel.mod

%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c -o $*.o $<

clean:
	rm -f dysurf validate_psf2d_parametrization *.o *.mod

.PHONY: clean
