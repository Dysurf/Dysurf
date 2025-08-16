from math import tanh, exp, ceil
import numpy as N
#from DebyeWallerCalculator import DebyeWallerCalculator
#from DebyeWallerCalculator import calcBoseEinstein
#from DW_cython import calcBoseEinstein
#from idf.Polarizations import read as readIDFpolarizations
#from idf.Omega2 import read as readIDFomega2s
import kernelGenerator2.SqeCalculator
from kernelGenerator2.SqeCalculator import *

class SQECalculatorForExistingPhononsForGivenQ:
    def __init__(self,SqeCalc) :
        #SQECalc is an instance of SQECalculator with associated values
        self._name='Class to calculate SQE for existed phonons given an SQE Object'
        self.SqeCalc=SqeCalc
        

    def calcSqeCohforExistedPhononsforGivenQ(self, Q):
            
        #sqe = 0
        SandE=[]
        qtransfer = Q
        #------------------------------
        # Calcualte Debye-Waller Factor
        #------------------------------
        d=0.0
        self.SqeCalc._DebyeWallerFactorList1 = []
        for atomIndex in range(len(self.SqeCalc._unitcell.getAtoms())):
            #print self.SqeCalc._displacements[atomIndex]
            dd= 0.0
            for j in range(3):
                #print self.SqeCalc._displacements[atomIndex][j+1]
                d= self.SqeCalc._displacements[atomIndex][j+1]*qtransfer[j]
                dd=dd+d
            dd=0.5*dd*dd
            self.SqeCalc._DebyeWallerFactorList1.append(dd)
            #print "self.SqeCalc._DebyeWallerFactorList1", self.SqeCalc._DebyeWallerFactorList1
            #print "qtransfer=",qtransfer
        #print self.SqeCalc._DebyeWallerFactorList
        #--------------------------------------
        # index the qtransfer
        #-------------------------------------
        
        
        
        #####################
        # 2013-05-06 JLN
        # Found that this is indexing the q transfer by linear searching though the kpoints.
        # Changed so that the kpoints are a list in the calcSqeCohforExistedPhononsforGivenQ object
        # and the searching is done by the list index rather than a linear search through the object
        # Speedup is x2.5
        tol=1.0e-8
        diff=self.SqeCalc._kpts - qtransfer
        idx=N.where((abs(diff) < tol).all(-1))
        kptIndex=idx[0][0]
        
        
        kvec = self.SqeCalc._kpts[kptIndex]
        
        for modeIndex in range(self.SqeCalc._D*self.SqeCalc._numatoms):
            energy = self.SqeCalc._energies[kptIndex][modeIndex]
            # only consider this mode if the phonon energy is close to E:
            #if N.abs(energy - E) > self.SqeCalc._etransferTol: continue
            # "throw away" imaginary modes:
            if energy < 0: continue
            # mode within tolerances, go on:
            vec = self.SqeCalc._polvecs[kptIndex][modeIndex]
            # Inner loop: we need to sum over all the atoms:
            innersum = 0
            for atom, atomIndex, pos in zip(self.SqeCalc._unitcell.getAtoms(),
                                            range(len(self.SqeCalc._unitcell.getAtoms())),
                                            self.SqeCalc._unitcell.getPositions()):
                pol = vec[atomIndex]
                qdote = N.dot(qtransfer,pol)
                #
                cpos = N.dot(self.SqeCalc._unitcell.getCellVectors(),pos)
                qdotd = N.dot(qtransfer,cpos)
#                qdotd = 0.0
                weight = N.exp(1j * qdotd - self.SqeCalc._DebyeWallerFactorList1[atomIndex])
                
                
                weight *= qdote
                weight *= atom.coh_sct_length / N.sqrt(atom.mass)
                innersum += weight
                
            sqe = (innersum.real**2 + innersum.imag**2) * (calcBoseEinstein(self.SqeCalc._temperature, energy) + 1.0) / energy
            SandE.append([energy,kvec,sqe])
        
        return SandE
        pass # end of calcSqeCoh

#        l=-1
#        for k in self.SqeCalc._kpts:
#            l=l+1
#            if k[0]==qtransfer[0] and k[1]==qtransfer[1] and k[2]==qtransfer[2]:
#               kptIndex= l
#               #print "l=", l
#               #print "qtransfer=",qtransfer
#               #print "k=", k 
#               #print "self.SqeCalc._kpts=", self.SqeCalc._kpts
#               #kptIndex = self.SqeCalc._kpts.index(qtransfer)
#               #print "kptIndex=", kptIndex
#               qred = self.SqeCalc.reduceQvector(qtransfer)
#               #print qred
#               kvec = self.SqeCalc._kpts[kptIndex]
#               #for kptIndex in range(len(self._kpts)):
#                   # check that the phonon wavevector is within neighborhood of qtransfer:
#                   # note: for the phonon annihilation process, we need to test that
#                   # qred is close to (-kvec)
#                   #kvec = self._kpts[kptIndex]
#                   #if N.dot((qred-kvec),(qred-kvec)) > self._qtransferTolRadius**2:
#                      #continue # skip to next phonon wavevector
#                      #print 'Found k-point in neighborhood of Q-transfer.'
#               for modeIndex in range(self.SqeCalc._D*self.SqeCalc._numatoms): 
#                   energy = self.SqeCalc._energies[kptIndex][modeIndex]
#                   # only consider this mode if the phonon energy is close to E:
#                   #if N.abs(energy - E) > self.SqeCalc._etransferTol: continue
#                   # "throw away" imaginary modes:
#                   if energy < 0: continue  
#                   # mode within tolerances, go on:
#                   vec = self.SqeCalc._polvecs[kptIndex][modeIndex]
#                   # Inner loop: we need to sum over all the atoms:
#                   innersum = 0
#                   for atom, atomIndex, pos in zip(self.SqeCalc._unitcell.getAtoms(),
#                                                range(len(self.SqeCalc._unitcell.getAtoms())),
#                                                self.SqeCalc._unitcell.getPositions()):
#                       pol = vec[atomIndex]
#                       qdote = N.dot(qtransfer,pol)
#                       #qdote = N.dot(wavevector,pol[:,0]) + np.dot(wavevector,pol[:,1]) * 1j
#                       # dot product of qtransfer with atom position in cell,
#                       # this requires that the qtransfer has been translated
#                       # into the reciprocal space of the crystal appropriately,
#                       # according to the orientation of the crystal
#                       qdotd = N.dot(qtransfer, pos)
#                       #weight = N.exp(1j * qdotd - self.SqeCalc._DebyeWallerFactorList1[atomIndex])
#                       #weight = N.exp(1j * qdotd ) # set debye waller factor to 0 
#                       weight = N.exp(1j - self.SqeCalc._DebyeWallerFactorList1[atomIndex]) #remove phase factor from position
#                       weight *= qdote
#                       weight *= atom.coh_sct_length / N.sqrt(atom.mass)
#                       innersum += weight
#               
#                   sqe = (innersum.real**2 + innersum.imag**2) * (calcBoseEinstein(self.SqeCalc._temperature, energy) + 1.0) / (energy) 
#                   SandE.append([energy,kvec,sqe])
#                               
#        return SandE
#        pass # end of calcSqeCoh

