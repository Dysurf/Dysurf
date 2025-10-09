.SUFFIXES: .f90

#-----------------------------------------------
# Makefile for the program 'Dysurf'.
# Please modify the variables properly.
# We recommend to use Intel fortran compiler.
#-----------------------------------------------

FC = ifx
FFLAGS = -traceback -debug -O2 -static-intel

FCL = ${FC}
LDFLAGS = -L/usr/lib64 -qmkl
#LDFLAGS= -llapack -lblas


LAPACK = 
LIBS = $(LAPACK)

#-----------------------------------------------
# General rules 
#-----------------------------------------------

FCSRC= $(wildcard *.f90)

MODS = ${FCSRC:.f90=.mod}

OBJS= ${FCSRC:.f90=.o}

default: dysurf

all: dysurf

dysurf: ${OBJS}
	${FCL} ${LDFLAGS} -o $@ $^ ${LIBS}
phonon_spectra.o phonon_spectra.mod: variables.mod constants.mod readin.mod func.mod
qpoints.o qpoints.mod: variables.mod
readin.o readin.mod: variables.mod constants.mod func.mod
rmsd.o rmsd.mod: phonon_spectra.mod func.mod constants.mod variables.mod
sqecalc.o sqecal.mod: func.mod variables.mod constants.mod
dysurf.o: constants.mod func.mod variables.mod phonon_spectra.mod qpoints.mod readin.mod rmsd.mod sqecalc.mod

%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod

clean:
	rm -f dysurf *.o *.mod

.PHONY: clean
