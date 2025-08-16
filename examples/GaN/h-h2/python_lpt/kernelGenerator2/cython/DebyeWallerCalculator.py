# DebyeWallerCalculator
# Olivier Delaire
# B. Keith
#from kernelGenerator.phononSqe.NetcdfPolarizationRead import NetcdfPolarizationRead

__doc__ = """Implementation of a Debye-Waller calculator.
Calculates the Debye-Waller factor for each atom in the unit cell,
based on a list of phonon modes passed as input."""

#from Units import *
import numpy as np
from DW_cython import getDWFactorForAtom_cython

#############################

class DebyeWallerCalculator:
    """Implementation of a Debye-Waller calculator.
    Calculates the Debye-Waller factor for each atom in the unit cell,
    based on a list of phonon modes passed as input."""

    hbar=6.58211899e-16*1e12*1e3#eV*s * ps/s * meV/eV = meV*ps

    def __init__(self, unitcell=None, kptlist=None, energies=None, polvecs=None, eigenVecFile=None):
        self._unitcell = unitcell
        self._kptlist = kptlist
        self._energies = energies
        self._polvecs = polvecs
        self.eigenVecFile=eigenVecFile
#        if len(polvecs) == 0:
#            raise ValueError, "DebyeWallerCalculator needs the phonon modes at least for one k-points."
        pass # end of __init__

    def getDWFactorAllAtoms(self, wavevector, temperature):
        """Returns the D-W factors for all the atoms in the unit cell at given wavevector transfer.
        The wavevector transfer is expected to be in same units as phonon wavevectors.
        The temperature is expected in Kelvin."""
        
        wavevector = np.array(wavevector)
#        T = temperature 
        DW = []
#        nkpt = len(self._kptlist)
#        # nphonons : number of points in the BZ is equal to number of uc's in crystal
#        kptindex = 0
#        natom = self._unitcell.getNumAtoms()
#        nModes = 3 * natom
#        weight = 0 # DW factor contribution
        atomindex = 0
        for atom in self._unitcell:
            # !!!
            # we have to make sure that the order of the atoms in the polarization vectors
            # is the same as the order in which the atoms are returned from the unit cell
            DW[atomindex] = self.getDWFactorForAtom(atomindex, wavevector, temperature)
            atomindex += 1
            pass
        # end of loop on atoms 

        return DW
    # enf of getDWFactorAllAtoms

    def getDWFactorForAtom(self, atomindex, wavevector, temperature):
        """Returns the D-W factor for an atom in the unit cell 
        (corresponding to atomindex),
        at given wavevector.
        Wavevector is expected to be in same units as phonon wavevectors.
        Temperature is expected in Kelvin."""
        wavevector = np.array(wavevector)
        
	DW = getDWFactorForAtom_cython(atomindex, self._unitcell[atomindex].getAtom().mass, len(self._kptlist), self._unitcell.getNumAtoms(), wavevector, self._energies, self._polvecs.real, self._polvecs.imag, temperature)
	
        return DW
      
