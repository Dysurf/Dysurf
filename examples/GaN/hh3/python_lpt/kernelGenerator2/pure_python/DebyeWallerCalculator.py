# DebyeWallerCalculator
# Olivier Delaire
# B. Keith
#from kernelGenerator.phononSqe.NetcdfPolarizationRead import NetcdfPolarizationRead

__doc__ = """Implementation of a Debye-Waller calculator. Pure Python code.
Calculates the Debye-Waller factor for each atom in the unit cell,
based on a list of phonon modes passed as input."""

#from Units import *
import numpy as np

def calcBoseEinstein(T,E,mu=0):
    """Helper function to calculate the Bose-Einstein distribution for energy E at temperature T.
    The chemical potential mu is set to zero by default. Pure Python code.
    Units are hard-coded as follows:
    Temperature: Kelvin
    Energy: meV
    mu: meV (chemical potential)
    """
    k_B = 8.617343E-2 # meV/K
    hbar = 6.582119E-13 # meV * s
    eps = 1e-6
    if (T<=0):
        raise ValueError, "Temperature must be positive."
    else:
        if (E == mu):
            # raise ValueError, "Bose Einstesin distribution infinite in zero."
            #rt = 1.0 / (np.exp((eps)/ (k_B*T)) - 1.0)
            rt = 1.0e9
        else:
            rt = 1.0 / (np.exp((E-mu)/ (k_B*T)) - 1.0)
        return rt
    

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
        """Returns the D-W factor for an atom in the unit cell (corresponding to atomindex),
        and at given wavevector.
        Wavevector is expected to be in same units as phonon wavevectors.
        Temperature is expected in Kelvin."""
        
        wavevector = np.array(wavevector)
        T = temperature # temperature is expected ot be in Kelvin
        DW = 0
        nkpt = len(self._kptlist)
        # nphonons : number of points in the BZ is equal to number of uc's in crystal
        kptindex = 0
        natom = self._unitcell.getNumAtoms()
        nModes = 3 * natom
        weight = 0 # DW factor contribution
        
        #print 'energies',self._energies 
        
        if self.eigenVecFile!=None:   
            #open file
            polRead=NetcdfPolarizationRead(self.eigenVecFile)
            #create DW     
            #print len(self._kptlist), nModes
            #import sys
            #sys.exit()
            for kptindex in range(len(self._kptlist)):
                for modeIndex in range(nModes): 
                    energy = self._energies[kptindex][modeIndex]
                    if energy<0: continue  #"throw away" imaginary modes
                    vec = polRead.readVec(kptindex,modeIndex)
                    
                    #print 'kptindex',kptindex
                    #print 'modeIndex',modeIndex
                    #print 'energy',energy
                    #print vec
                    #print vec.shape
                    #print 'atom index',atomindex
                    pol = vec[atomindex]
                    #print pol
                    #weight = np.dot(wavevector,pol[:,0])**2 +  np.dot(wavevector,pol[:,1])**2
                    weight = np.dot(wavevector,pol.real)**2 +  np.dot(wavevector,pol.imag)**2
                    # weight is square-modulus of (real) wavevector dotted with complex polarization
                    thermalfactor = 2.0 * calcBoseEinstein(T, energy) + 1.0
                    #print 'weight',weight
                    #print 'thermalfactor', thermalfactor
                    if ( energy > 1e-6):
                        DW += weight * thermalfactor / energy
        else:
            for kptindex in range(len(self._kptlist)):
                for modeIndex in range(nModes): 
                    energy = self._energies[kptindex][modeIndex]
                    pol = self._polvecs[kptindex][modeIndex][atomindex]
                    #weight = np.dot(wavevector,pol[:,0])**2 +  np.dot(wavevector,pol[:,1])**2
                    weight = np.dot(wavevector,pol.real)**2 +  np.dot(wavevector,pol.imag)**2
                    # weight is square-modulus of (real) wavevector dotted with complex polarization
                    thermalfactor = 2.0 * calcBoseEinstein(T, energy) + 1.0
                    if ( energy > 1e-6):
                        DW += weight * thermalfactor / energy        
          
        mass = self._unitcell[atomindex].getAtom().mass
        # normalization, cf Squires (3.74)
        # (we normalize by nkpt, since in the number of kpoints in the BZ is equal
        # to the number of unit cells in the crystal)
        #print 'mass',mass
        #print 'hbar2nd',self.hbar
        meVps_to_uAng2Byps2=9.6485341
        DW *= (self.hbar*meVps_to_uAng2Byps2 / (4.0 * mass * nkpt) )
        return DW

def test():
    from UnitCell import UnitCell,Site
    from Atom import Atom

    uc = UnitCell()
    at0 = Atom(symbol="Al")
    at1 = Atom(symbol="Fe")
    site0=Site([0.0, 0.0, 0.0],at0)
    site1=Site([0.5, 0.5, 0.5],at1)
    uc.addSite(site0,"")
    uc.addSite(site1, "")
    print uc
    print

    from python.parsing.parse_phon_results import parse
    phonlist = parse("phon.out_dos_noIBZ")

    
    wvectors = [[qx, 0.0, 0.0] for qx in range(30)]

    dwc = DebyeWallerCalculator(uc, phonlist, phonlist)

    dwlist = [dwc.getDWFactorForAtom(0, wavevector, 300) for wavevector in wvectors]
    #print dwlist
    return dwlist
    
if __name__ == "__main__": test()
