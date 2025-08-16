#!/usr/bin/python
import pylab as pl
from crystal.UnitCell import UnitCell,Site
from crystal.SimpleAtom import Atom
from idf.Polarizations import read as readIDFpolarizations
from idf.Omega2 import read as readIDFomega2s
from DebyeWallerCalculator import DebyeWallerCalculator

def test():
    uc = UnitCell()
    at0 = Atom(symbol="Al", mass=26.982)
    at1 = Atom(symbol="Fe", mass=55.85)
    site0=Site([0.0, 0.0, 0.0],at0)
    site1=Site([0.5, 0.5, 0.5],at1)
    uc.addSite(site0,"")
    uc.addSite(site1, "")
    print uc
    poldata = readIDFpolarizations('FeAl_polarizations.idf')
    pols = poldata[1]
    edata = readIDFomega2s('FeAl_omega2s.idf')
    energies = edata[1]
    wvectors = [[qx, 0.0, 0.0] for qx in range(30)]
    dwc = DebyeWallerCalculator(uc, pols, energies, pols)
    dwlist = [dwc.getDWFactorForAtom(0, wavevector, 300) for wavevector in wvectors]
    pl.plot(dwlist)
    
if __name__ == "__main__": 
    test()
