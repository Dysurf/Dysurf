.SUFFIXES: .f90

#-----------------------------------------------
# Makefile for the program 'Dysurf'.
# Please modify the variables properly.
# We recommend to use Intel fortran compiler.
#-----------------------------------------------
# 编译器设置
USE_MPI ?= 0
ifeq ($(USE_MPI),1)
FC = /usr/bin/mpif90
else
FC = gfortran
endif
FCL = ${FC}

# 编译选项
FFLAGS = -g -O2 -ffree-line-length-none -fPIE

# 链接选项
LDFLAGS = -L/usr/lib/x86_64-linux-gnu/ -llapack -lopenblas -lgfortran -lm -lpthread

#-----------------------------------------------
# General rules 
#-----------------------------------------------

PROGRAM_SRCS = dysurf.f90 validate_psf2d_parametrization.f90
FCSRC = $(filter-out $(PROGRAM_SRCS),$(wildcard *.f90))

MODS = ${FCSRC:.f90=.mod}

OBJS = ${FCSRC:.f90=.o}
PROGRAM_OBJS = dysurf.o validate_psf2d_parametrization.o

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
