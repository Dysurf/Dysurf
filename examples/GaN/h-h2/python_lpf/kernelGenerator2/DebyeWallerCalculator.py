# DebyeWallerCalculator
# Olivier Delaire
# B. Keith
#from kernelGenerator.phononSqe.NetcdfPolarizationRead import NetcdfPolarizationRead

__doc__ = """Implementation of a Debye-Waller calculator.
Calculates the Debye-Waller factor for each atom in the unit cell,
based on a list of phonon modes passed as input."""

#from Units import *
import numpy as np
import math
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

    def getDWFactorForCubicSymmetry(self,temperature):
        atomindex=0
        hbarsquare= 1.112121e-68 # J^2.s^2
        mass=self._unitcell[atomindex].getAtom().mass
        print ("mass=", mass)
        #mass = 92.906
        #T=T   # J
        M=1.6605389e-27  #kg
        q=1e+10
        cons= hbarsquare*q**2/(6*M)
        #cons= hbarsquare/(6*M)
        sumD=0
        DW=0
        k=0
        for i in range(48400):
            for j in range(3):
            #energies[i][j]=energies[i][j]  # ene in J
            #if self._energies[i][j]<0: continue
            #if self._energies[i][j]>0:
                k=k+1
                coth= 1/(math.tanh(0.5*self._energies[i][j]*11.594/temperature))
                #print energies[i][j], coth
                #sumD +=coth*6.241509752e+21/self._energies[i][j]
                sumD +=coth*1.60217646e-22/self._energies[i][j]
        DW=sumD*cons
        print ("cons=", cons)
        
        #DW=sumD*cons/k
        print ("DWcubic=", DW)
        self._DWD = DW
        return self._DWD


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
        print ("DW-non=", DW)
        return DW
      
