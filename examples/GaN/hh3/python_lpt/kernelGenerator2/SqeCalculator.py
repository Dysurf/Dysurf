#from crystal import CrystalStructure
#2
from math import tanh, exp, ceil
import numpy as N
from kernelGenerator2.DebyeWallerCalculator import DebyeWallerCalculator
#from DebyeWallerCalculator import calcBoseEinstein
from DW_cython import calcBoseEinstein
from idf.Polarizations import read as readIDFpolarizations
from idf.Omega2 import read as readIDFomega2s

class SqeCalculator:
    """A class that calculate the neutron scattering function,  S(q_vector, E),
    for scattering by phonons.
    The computation of S(q_vector, E) uses phonon eigenvectors (polarizations)
    and eigenvalues (phonon energies)."""
    
    gammaPointCalc = False
    hbar = 6.58211899e-16*1e12*1e3#eV*s * ps/s * meV/eV = meV*ps
    kB = 8.61734315e-5*1e3 #ev/K * meV/eV
    hbar_meVs = 6.582119*1e-13 # meV * s
    
    #eventually this can be generalized to liquids or polymers
    def __init__(self,
                 unitcell,
                 phononEnergies=None,
                 phononFrequencies=None,
                 polarizations=None,
                 kpoints=None,
                 tau = None,
                 energies=None,
                 displacements=None,
                 temperature=300.0,
                 frequencyLimits=(0,20000),
                 eigenVecFile=None,
                 D=3):
        """
        unitcell: a UnitCell object describing the crystal structure
        energies: phonon energies
        polarizations: phonon polarization vectors
        kpoints: list of phonon k-points, in reciprocal lattice fractional coords
        tau: a reciprocal lattice vector (in Miller notation)
        """
        # material properties:
        self.phononFrequencies=phononFrequencies
        self._energies = energies
        self._polvecs = polarizations
        self._displacements = displacements
        self._numkpts = None
        self._kpts = kpoints
        if kpoints is None:
            self._numkpts = None
        else:
            self._numkpts = len(kpoints)
        #self.setKpoints(kpoints)
        self._tau = None
        self._D = D
        self._unitcell = unitcell
        self._numatoms = None      # this gets set to the number of atoms in the cell next line: suppress this line?
        self.setUnitCell(unitcell)
        self.eigenVecFile=eigenVecFile
        #self.frequencyLimits=frequencyLimits

        # scattering properties:
        self._temperature = temperature
        self._etransfer = None
        self._qtransfer = None
        self._DebyeWaller = None
        self._DebyeWallerFactorList = None
        self._DebyeWallerFactorList1 = None
        # various properties, computation tolerances, etc.
        self._qtransferTolRadius = 0.1 # small number representing radius of neighborhood tolerance around qtransfer
        self._etransferTol = 1.0  # tolerance on the energy transfer (in meV)
        self.binSize=1.0

        # put all the eigenvalues into bins
        # count the bins from zero and go by binsize
        # but plot the bins vs the bin centers

        frequencyRange=int(ceil((frequencyLimits[1]-frequencyLimits[0])/self.binSize))
        ##self.binCounts=[0 for i in range(frequencyRange)]
        #self.binnedEigs=[[] for i in range(frequencyRange)]
        #self.binCenters=[i+self.binSize/2.0 for i in range(frequencyRange)]
        #for kptIndex in range(len(self._kpts)):
        #    for modeIndex in range(len(self.phononFrequencies[0])):
        #        if self.phononFrequencies[kptIndex,modeIndex] < 0: continue # ignore negative frequencies
                #int (eig/binSize) gives bin number
        #        binNumber=int(self.phononFrequencies[kptIndex,modeIndex]/self.binSize)
                #print binNumber
        #        self.binnedEigs[binNumber].append((self.phononFrequencies[kptIndex,modeIndex], kptIndex, modeIndex))
            
#        self.sIndex=zip(self._kpts,range(len(self._energies)))

#        for eig in self.binnedEigs:
#            print 'eig ', eig


        #self.polRead=NetcdfPolarizationRead(self.eigenVecFile)
        
        #self.dwc = DebyeWallerCalculator(self._unitcell, self._kpts, self.phononFrequencies, 
        #   eigenVecFile=self.eigenVecFile)

        self._DebyeWallerCalculator = DebyeWallerCalculator(self._unitcell,
                                                            self._kpts,
                                                            self._energies,
                                                            self._polvecs)
                                                            
    def calcSqeCohCreateAllmodes(self, Q, E):
        """Calculates the coherent S(q_vector,E), from the phonon eigenvalues and eigenvectors,
        at the q-pts stored and for the energies defined in self._energies.
        Based on Squires, Eq. (3.120)"""

        sqe = 0
        qtransfer = Q
        # first we calculate the DW factor for all the atoms in the unit cell at wavevector Q:
        self._DebyeWallerFactorList = []
        for atom, atomIndex in zip(self._unitcell.getAtoms(),range(len(self._unitcell.getAtoms()))):
            DWd=self._DebyeWallerCalculator.getDWFactorForAtom(atomIndex, Q, self._temperature)
            self._DebyeWallerFactorList.append(DWd)

        # Outer loop on all the modes: weed out the modes not within neighborhood of qtransfer
        # First, we reduce the Qtransfer vector into Brillouin Zone
        qred = self.reduceQvector(qtransfer)

        for kptIndex in range(len(self._kpts)):
            # check that the phonon wavevector is within neighborhood of qtransfer:
            # note: for the phonon annihilation process, we need to test that
            # qred is close to (-kvec)
            kvec = self._kpts[kptIndex]
            if N.dot((qred-kvec),(qred-kvec)) > self._qtransferTolRadius**2:
                continue # skip to next phonon wavevector
            #print 'Found k-point in neighborhood of Q-transfer.'
            for modeIndex in range(self._D*self._numatoms): 
                energy = self._energies[kptIndex][modeIndex]
                # only consider this mode if the phonon energy is close to E:
                if N.abs(energy - E) > self._etransferTol: continue
                # "throw away" imaginary modes:
                if energy < 0: continue  
                # mode within tolerances, go on:
                vec = self._polvecs[kptIndex][modeIndex]
                # Inner loop: we need to sum over all the atoms:
                innersum = 0
                for atom, atomIndex, pos in zip(self._unitcell.getAtoms(),
                                                range(len(self._unitcell.getAtoms())),
                                                self._unitcell.getPositions()):
                    pol = vec[atomIndex]
                    qdote = N.dot(qtransfer,pol)
                    #qdote = N.dot(wavevector,pol[:,0]) + np.dot(wavevector,pol[:,1]) * 1j
                    # dot product of qtransfer with atom position in cell,
                    # this requires that the qtransfer has been translated
                    # into the reciprocal space of the crystal appropriately,
                    # according to the orientation of the crystal
                    qdotd = N.dot(qtransfer, pos)
                    weight = N.exp(1j * qdotd - self._DebyeWallerFactorList[atomIndex])
                    weight *= qdote
                    weight *= atom.average_neutron_coh_xs / N.sqrt(atom.mass)
                    innersum += weight

                sqe += (innersum.real**2 + innersum.imag**2) * (calcBoseEinstein(self._temperature, energy) + 1.0) / energy
                    
        return sqe
        pass # end of calcSqeCoh
    
    
        
    def calcSqomegaIncoherent(self,Q,omega):
        '''This comes from Squires, combining equations 3.138 and a modification of 4.13.
        
Note that we have replaced the delta function with a gaussian.        
Since this goes to zero if omega should be in ps.'''

        negligibleFactor=10**-12

#        wvectors = [[qx, 0.0, 0.0] for qx in range(10)]
#
#        dwlist = [self.dwc.getDWFactorForAtom(0, wavevector, 300) for wavevector in wvectors]
                
        # sum over the  atoms in the unit cell
        #atomIndex=0
        sqomegaSum=0.0
        for atom, atomIndex in zip(self._unitcell.getAtoms(),range(len(self._unitcell.getAtoms()))):
            #loop over the frequencies that are within the same bin as the frequency given in the function arguments (omega)
            omegaBin=int(ceil(omega/self.binSize))
            deltaFunctionEigs=self.binnedEigs[omegaBin]
            innerSum=0.0
            for omegaPacket in deltaFunctionEigs:
                # unpacking
                omega,kptIndex,modeIndex = omegaPacket
                vec = self.polRead.readVec(kptIndex,modeIndex)
                polarizationVector = vec[atomIndex]
                #print pol
                modulusSquared = N.dot(Q,polarizationVector[:,0])**2 +  N.dot(Q,polarizationVector[:,1])**2

                innerSum+=modulusSquared/omega*self.boseFactorPlus1(omega)
                #print 'modulusSquared', modulusSquared
                #print 'innerSum', innerSum
                
                if innerSum > negligibleFactor:
                    print ('modulusSquared', modulusSquared)
                    print ('innerSum', innerSum)
            
            if innerSum < negligibleFactor: continue
            #mass=atom.mass
            #print self.dwc.getDWFactorForAtom(atomIndex, Q, self._temperature)
            DWd=self._DebyeWallerCalculator.getDWFactorForAtom(atomIndex, Q, self._temperature)
#            print atom.mass
#            print atom.average_neutron_inc_xs
            prefactor=1/atom.mass*atom.average_neutron_inc_xs**2*exp(-2.0*DWd)
            
            sqomegaSum+=prefactor*innerSum
            #atomIndex+=1
            
            
        return sqomegaSum
        
    def beta(self):
        return 1./(self.kB*self._temperature)
        
    def boseFactorPlus1(self,omega):
        trig=tanh(0.5*self.hbar*omega*self.beta())
        return 0.5*((1+trig)/trig)

    def setUnitCell(self, unitcell):
        """Sets the unit cell."""
        try:
            self._numatoms = unitcell.getNumAtoms()
        except:
            raise ValueError( 'unitcell should be a UnitCell object.')

    def getUnitCell(self):
        """Returns the unit cell"""
        return self._unitcell
    
    def setKpoints(self, kpoints):
        """Sets the phonon k-points.
        These must correspond to the points at which the polarizations and energies have been calculated."""
        self._kpts = kpoints
        self._DebyeWallerCalculator._kptlist = kpoints
        self._numkpts = len(kpoints)

    def getPhononEnergies(self):
        """Returns the phonon energies"""
        return self._energies

    def setPhononEnergies(self, energies):
        """Sets the phonon eigenvalues"""
        self._energies = energies
        self._DebyeWallerCalculator._energies = energies

#    def getPolarizationVectors(self):
#        """Returns the phonon eigenvectors"""
#        return self.

    def setPolarizationVectors(self, polvecs):
        """Sets the phonon eigenvectors"""
        self._polvecs = polvecs
        self._DebyeWallerCalculator._polvecs = polvecs

    def setTau(self, tau):
        """sets the reciprocal lattice vector tau,
        corresponding of to a shift of phonon q-points."""
        self._tau = tau

    def getTau(self):
        """Returns the reciprocal lattice vector tau,
        corresponding of to a shift of phonon q-points."""
        return self._tau
    
    # reading/writing of IDF data...    

    def readEigenvaluesFromIDFomega2(self, filename='omega2.idf'):
        """Reads in the phonon eigenvalues from an IDF-format file."""
        #try:
        idfdata = readIDFomega2s(filename=filename)
        #except:
        #    raise IOError, 'Could not read the IDF data for Omega2.'
        infotuple = idfdata[0]
        print (infotuple)
        om2 = idfdata[1]
        # check the dimensions of the data:
        print (om2.shape)
        if om2.shape != (self._numkpts, self._numatoms*self._D):
            print ('The eigenvalues from file '+filename+' do not have the proper dimensions.')
            print ('Check number of q-points and number of atoms.')
            return
        energies = N.zeros((self._numkpts, self._numatoms*self._D), dtype='float')
        # We now convert the angular frequencies squared into energies (meV)
        # E_meV = THz2meV * sqrt(om2_THz^2)
        THz2meV = 4.1357 # meV * s
        for kIndex in range(self._numkpts):
            for modeIndex in range(self._numatoms * self._D):
                energies[kIndex][modeIndex] = THz2meV * om2[kIndex, modeIndex]
        self.setPhononEnergies(energies)
        return

    def writeIDFomega2FromEigenvalues(self, filename='omega2.idf'):
        """Writes the phonon eigenvalues to an IDF-format file."""
        om2 = N.zeros((self._numkpts, self._numatoms*self._D), dtype='float')
        # We now convert the energies (meV) into angular frequencies squared (Hz^2)
        # om2 = ( E(meV) / hbar_meVs )^2
        hbar_meVs = 6.582119*1e-13 # meV * s
        for kIndex in range(self._numkpts):
            for modeIndex in range(self._numatoms * self._D):
                om2 = ( self._energies[kIndex, modeIndex] / hbar_meVs)**2
        idf.Omega2.write(om2, filename=filename, comment=str(self._unitcell), D=self._D)
        return

    def readIDFeigenvectors(self, filename='polarizations.idf'):
        """Reads in the phonon eigenvectors from an IDF-format file."""
        #try:
        idffilename = filename
        idfdata = readIDFpolarizations(filename=idffilename)
        #except:
        #    raise IOError, 'Could not read the IDF data for Polarizations.'
        #idffilename = filename
        #idfdata = idf.Polarizations.Read(filename=idffilename)
        infotuple = idfdata[0]
        print (infotuple)
        self.setPolarizationVectors(idfdata[1])
        
    def writeIDFeigenvectors(self, filename='polarizations.idf'):
        """Writes the phonon eigenvectors to an IDF-format file."""
        try:
            idf.Polarizations.write(self._polvecs,
                                    filename=filename,
                                    comment=str(self._unitcell))
        except:
            raise IOError ('Could not write the Polarizations to IDF-format file.')
        

    def calcSqeInc(self):
        """Calculates the incoherent S(q_vector,E), from the phonon eigenvalues and eigenvectors,
        at the q-pts stored and for the energies defined in self._energies"""
        raise NotImplementedError


    def calcDebyeWaller(self):
        """Calculates the Debye-Waller factor W_d for all the atoms in the unit cell."""
        raise NotImplementedError

    def calcDebyeWallerForAtom(self, atomidex):
        """Calculate the Debye-Waller factor W_d for all one atom."""
        raise NotImplementedError
    
    def setBinSize(self,bin):
        self.bin=bin


    def reduceQvector(self, q):
        """Brings cartesian reciprocal vector q=(x,y,z) 'into' the reciprocal cell in which phonons are computed,i.e.: (x,y,z)->(x',y',z') such that x',y',z' in [0,1( """
        uc = self._unitcell
        ucstar = uc.getReciprocalUnitCell()
        ucstar.setCellVectors(ucstar.getCellVectors() / 2. / N.pi)
        qfrac = ucstar.cartesianToFractional(q)
        qred = ucstar.bringFractionalPositionIntoCenteredCell(qfrac)
        qredcart = ucstar.fractionalToCartesian(qred)
	#qred = N.array(q)
        #assert (len(qred) == 3)
        #for i in range(3):
        #    if qred[i]<0:
        #        while qred[i]<0:
        #            qred[i] += 1
        #    if qred[i]>=1:
        #        while qred[i]>=1:
        #            qred[i] -= 1
        return qredcart



##     def calcIntensityGrid(self):
##         """This is a helper function to calculate the scattering intensity,
##         from the polarization vectors and the scattering vector Q = q + tau.
##         It requires that the polarization vectors have been set,
##         for a given atom and branchindex.
##         The reciprocal lattice vector tau is a shift in reciprocal space."""
##         if self._tau is not None:
##             qshift = self._tau
##         else:
##             qshift = [0,0,0]

##         cellvectors = self._unitcell.getCellVectors()

##         space = VectorSpaces.VectorSpaceWithBasis(cellvectors.tolist())
##         origin = Vector(np.array([0,0,0]))
##         intgrid = Grid(space=space, origin=origin)

##         #kptarray = self._kptGrid.GetArray()
##         #polarray = self._polarizationGrid.GetArray()

##         # nkpt = len(self._phonondata)
##         intensity = N.zeros(self._numkpts, dtype=float)

##         for qindex in range(self._numkpts):
##             # for a non-monatomic crystal we only calculate the orientation-dependent part of the
##             # scattering intensity for each atom and branch.
##             # The relative intensities for scattering from different atoms are not preserved.
##             # We need to compute the Debye-Waller factor for each atom to get this.
##             # Here we do the dot-product with the "full" imaginary polarization vector,
##             # since the real polarization vector is obtained by a unitary transfo, wich
##             # does not change the value of | q_vec . e_vec |^2
            
##             dot = N.dot(( kshift + self._qpts[qindex]), self.+polvecs[qindex])
##             z = dot[0] + dot[1] * 1j
##             intensity[k0][k1][k2] = abs(z)**2

##         intgrid.SetArray(intensity)
##         self._intensityGrid = intgrid

##         return


