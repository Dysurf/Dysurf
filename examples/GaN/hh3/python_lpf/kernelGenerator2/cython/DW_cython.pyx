cimport numpy as np
ctypedef np.float64_t DOUBLE
import numpy

#############################

def calcBoseEinstein(float T,float E,float mu=0):
    """Helper function to calculate the Bose-Einstein distribution for energy E at temperature T.
    The chemical potential mu is set to zero by default.
    Units are hard-coded as follows:
    Temperature: Kelvin
    Energy: meV
    mu: meV (chemical potential)
    """
    cdef float k_B = 8.617343E-2 # meV/K
    cdef float hbar = 6.582119E-13 # meV * s
    cdef float eps = 0.0 #1e-6
    cdef float rt
    if (E == mu):
        rt = 1.0 / (numpy.exp((eps)/ (k_B*T)) - 1.0)
    else:
        rt = 1.0 / (numpy.exp((E-mu)/ (k_B*T)) - 1.0)
    return rt

#############################

def getDWFactorForAtom_cython(int atomindex, float atommass, int nkpt, int natom,  np.ndarray[DOUBLE, ndim=1] wavevector, np.ndarray[DOUBLE, ndim=2] energies, np.ndarray[DOUBLE, ndim=4] polvecs_real, np.ndarray[DOUBLE, ndim=4] polvecs_imag, float temperature):
    """Returns the D-W factor for an atom in the unit cell (corresponding to atomindex),
    and at given wavevector.
    Wavevector is expected to be in same units as phonon wavevectors.
    Temperature is expected in Kelvin."""

    cdef float T = temperature
    cdef float DW = 0.0
    cdef float hbar = 6.582119E-13
    cdef int kptindex = 0
    cdef int modeindex = 0
    cdef int nModes = 3 * natom
    cdef float weight = 0.0
    cdef float thermalfactor = 0.0
    cdef float energy = 0.0
    #cdef np.ndarray[DOUBLE] pol_real
    #cdef np.ndarray[DOUBLE] pol_imag

    #for kptindex in range(nkpt):
        #for modeIndex in range(nModes):
            #energy = energies[kptindex][modeIndex]
            #pol_real = polvecs_real[kptindex][modeIndex][atomindex]
            #pol_imag = polvecs_imag[kptindex][modeIndex][atomindex]
            #weight = numpy.dot(wavevector,pol_real)**2 +  numpy.dot(wavevector,pol_imag)**2
            # weight is square-modulus of (real) wavevector dotted with complex polarization
            #thermalfactor = 2.0 * calcBoseEinstein(T, energy) + 1.0
            #DW += weight * thermalfactor / energy

    weights = numpy.dot(polvecs_real[:,:,atomindex,:], wavevector)**2 + numpy.dot(polvecs_imag[:,:,atomindex,:], wavevector)**2
    thermalfactors = 1.0 / (numpy.exp( energies / (8.617343E-2 * T)) - 1.0)
    DW = (weights * thermalfactors / energies).sum()

    # normalization, cf Squires (3.74)
    # (we normalize by nkpt, since in the number of kpoints in the BZ is equal
    # to the number of unit cells in the crystal)
    cdef meVps_to_uAng2Byps2 = 9.6485341
    DW *= (hbar*meVps_to_uAng2Byps2 / (4.0 * atommass * nkpt) )
    return DW

