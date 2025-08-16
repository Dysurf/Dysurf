import numpy as np
import numpy.linalg as la
from crystal.SimpleAtom import Atom
from crystal.crystalUtils.MonkhorstPack import MonkhorstPack


##########################################################

class Site:
    """Representation of a crystallographic site.
    A site corresponds to a position in fractional coordinates.
    It also may have an atom associated with it."""

    def __init__(self, position=(0.0,0.0,0.0), atom=None, occproba=1.0):
        for x in position:
            if (x < 0) or (x > 1):
                raise ValueError ("Site coordinates must be fractional positions.")
        self._position = np.array( position )
        self._atom = atom
        self._occproba = occproba
        self.xyz=self._position

    def __str__(self):
        rt = str(self._position) + ":" + str(self._atom)
        return rt

    def setPosition(self, position):
        """Sets the position (in fractional coordinates) of a site.""" 
        for x in position:
            if (x < 0) or (x > 1):
                raise ValueError("Site coordinates must be fractional positions.")
        self._position = np.array(position)
        
    def getPosition(self):
        """Get the position (in fractional coordinates) of a site."""
        return np.array(self._position)

    def setAtom(self, atom):
        """Set an atom at a site."""
        self._atom = atom

    def getAtom(self):
        """Returns the atom at a site."""
        return self._atom

    def getOccProba(self):
        """Returns the occupation probability for the atom at this site."""
        return self._occproba

    def setOccProba(self,p):
        """Sets the occupation probability for the atom at this site."""
        try:
            proba = float(p)
        except:
            raise ValueError( 'Probability should be a number in [0,1].')
        if proba <= 1.0 and proba >= 0.0:
            self.__occproba = proba
        else:
            raise ValueError ('Probability should be a number in [0,1].')
            
        
    pass # end of class site

##########################################################

class UnitCell:  
    """Representation of a crystal unit cell."""
#    unitCellWrapper=uniCellWrapper

# We do not want to make UnitCell a component,
# since it represents data, rather than a process.
        #atoms = multiLineStr('Atoms',default='''"H" 0.0 0.0 0.0\n"H" 1.0 0.0 0.0''')  
        #atomFile = inv.str('atom and position file',defa)
        #atomFile
#        for i in range(UnitCell.numAtoms):
#            line1 = '''atom%s=pinv.list( atom%s, default = "'H' 0.0 0.0 0.0")''' % (i,i)
#            line2 = "atom%s.meta['tip'] = 'species and position of atom %s'" % (i,i)
#            exec line1 in locals()
#            exec line2 in locals()    
    
    def __init__(self, name='UnitCell',cellvectors=None, spaceGroup=None):
        if cellvectors is None:
            cellvectors = np.array( [ (1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0) ] )

        self._cellvectors = np.array(cellvectors)
        self._spaceGroup = spaceGroup
        self._sites = []  # list of sites
        self._siteIds = {}  # dictionary {siteId : site }   
        self._n = 0
        self._numAtoms = 0
        self._numSites = 0
        self.base=self._cellvectors
        return

    def __iter__(self): return self._sites.__iter__()

    def __getitem__(self, i): return self._sites[i]

    def __str__(self):
        rt = "UnitCell a=%s, b=%s, c=%s\n" % tuple(self._cellvectors)
        for siteId in self._siteIds.keys():
            rt += "\n%s\n" % siteId
            rt += "\n position: %s\n" % (self._siteIds[siteId].getPosition()) 
            rt += "\n%s \n" % (self._siteIds[siteId].getAtom())
            continue
        return rt
    
    def __len__(self):
        return self.getNumAtoms()

    def __copy__(self):
        new = UnitCell()
        new._cellvectors = self._cellvectors
        for siteId in self._siteIds.keys():
            newsite = Site(self._siteIds[siteId].getPosition(),
                           Atom(Z=self._siteIds[siteId].getAtom().Z))
            new.addSite(newsite, siteId = siteId)
        return new


    def addSite(self, site, siteId=None):
        """Adds a site to the unit cell."""
        #assert ( isinstance(site, Site) )

        if siteId.__class__ is not "string".__class__ :
            raise ValueError( 'site Id should be a string!')
        pass
        
        if siteId is None or siteId == "": siteId= "Id%s" % (self._uniqueID() )

        self._sites.append( site )
        self._siteIds[siteId] = self._sites[-1]
        self._numSites += 1
        if site.getAtom() is not None:
            self._numAtoms +=1
        return

    def addAtom(self, atom, position, siteId=None):
        """add( Fe_atom,  np.array((0.25,0,0)) ) --> adds atom Fe_atom to UnitCell
        add(Atom(Z=26), (0.25, 0, 0)) --> adds atom Fe_atom to UnitCell
        """
        assert ( isinstance(atom, Atom) )

        if siteId.__class__ is not "string".__class__ :
            raise ValueError ( 'site Id should be a string!')
        pass 
        
        if siteId is None or siteId == "": siteId= "Id%s" % (self._uniqueID() )

        newSite = Site(position, atom)
        self._sites.append( newSite )
        self._siteIds[siteId] = self._sites[-1]
        self._numAtoms += 1
        return

    def getNumAtoms(self):
        """Returns the number of atoms in the unit cell."""
        return self._numAtoms

    def getNumSites(self):
        """Returns the number of sites in the unit cell."""
        return self._numSites

    def getAtomTypeDenum(self):
        """Returns a dictionary with the types of atoms (symbols) and the number of atoms of each type."""
        denum = {}
        symbolmass = [(s.getAtom().symbol, s.getAtom().mass) for s in self._sites]
        for s in symbolmass:
            if not denum.has_key(s):
                denum[s] = symbolmass.count(s)
        return denum

    def getAtoms(self):
        """Returns a list of atoms for all sites in the unit cell.
        Sites with no atom return none."""
        return [s.getAtom() for s in self._sites]

    def getPositions(self):
        """Returns a list of positions for all sites in the unit cell."""
        return [s.getPosition() for s in self._sites]

    def getProperties(self, siteId):
        """ Returns the properties of the atom at a site."""
        return str( self._siteIds[siteId].getAtom() )

    def getPosition(self, siteId):
        """Returns the (fractional) position of a site."""
        return self._siteIds[siteId].getPosition()

    def setPositions(self, positions):
        """Sets the (fractional) positions of the sites in the unit cell."""
        assert(len(positions) == self.getNumSites())
        for isite in range(self.getNumSites()):
            self._sites[isite].setPosition(positions[isite])

    def getCartesianPosition(self, siteId):
        """Returns the cartesian position of a site."""
        return self.fractionalToCartesian(self._siteIds[siteId].getPosition())

    def fractionalToCartesian(self, fracpos):
        """Converts a coordinate from fractional to cartesian.""" 
        return (fracpos * self._cellvectors).sum(1)  # should be double-checked

    def cartesianToFractional(self, cartpos):
        """Converts a coordinate from cartesian to fractional."""        
        """ should be double-checked"""
        return (cartpos * la.inv(self._cellvectors)).sum(1) 
    
    def getIds(self):
        """Returns the list of Ids for sites in the unit cell."""
        return self._siteIds.keys()

    def getSite(self, siteNum):
        """Returns the Site number 'num'."""
        return self.__getitem__(siteNum)

    def getSiteFromId(self, siteId):
        """Returns a Site instance from its Id in the unit cell."""
        return self._siteIds[siteId]


    def _uniqueID(self):
        rt = self._n
        self._n += 1
        return rt

    def getCellVectors(self):
        """Returns the UnitCell cell vectors."""
        return self._cellvectors

    def setCellVectors(self, cellvecs):
        """Set the unit cell vectors (lattice vectors)."""
        # should add some checking here
        self._cellvectors = np.array(cellvecs)
        return

    def getVolume(self):
        """
        Returns the volume of the unit cell: |det(a1, a2, a3)|.
        Uses Numpy.linalg."""
        return abs(la.det(self._cellvectors))

    def getRecipVectors(self):
        """Returns the reciprocal lattice vectors (with the 2pi factors)."""
        recipvectors = 2 * np.pi * la.inv(np.transpose(self._cellvectors))
        # this needs to be checked
        return recipvectors

    def getReciprocalUnitCell(self):
        """
        Returns the reciprocal space unit cell, in the
        crystallographic sense (i.e. with the 2pi factors).
        Uses Numpy.linalg."""
        
        recipvectors = 2 * np.pi * la.inv(np.transpose(self._cellvectors))
        recipUC = UnitCell(cellvectors=recipvectors)
        return recipUC

    def computeDistances(self, maxdist=30, latticeRange=[2,2,2]):
        """ unitcell.computeDistances(self, [nx,ny,nz]):
        builds up a Big multiple dictionary, namely
        self.distances[atom1][atom2][(DX,DY,DZ)]
        (DX,DY,DZ) are integer numbers specifying the cell containing atom2,
        relatively to atom1.
        DX,DY,DZ run from -nx to nx, -ny to ny, -nz to nz, respectively."""
        distances = {}
        idlist = self.getIds()
        for iA in range(0, len(idlist)):
            idA = idlist[iA]
            distances[idA] = {}
            for iB in range(0, len(idlist)):
                idB = idlist[iB]
                distances[idA][idB]={}
                for tx in range(-latticeRange[0],latticeRange[0]+1):
                    for ty in range(-latticeRange[1],latticeRange[1]+1):
                        for tz in range(-latticeRange[2],latticeRange[2]+1):
                            posA = self.getCartesianPosition(idA)
                            posB = self.getCartesianPosition(idB) + np.dot([tx,ty,tz], self._cellvectors)
                            dist = np.sqrt(np.sum( (posB-posA) * (posB-posA) ))
                            if(dist<maxdist):
                                distances[idA][idB][(tx,ty,tz)] = dist
        self._distances = distances
        return

        
    def bringFractionalPositionIntoCell(self, fracpos):
        """Brings a fractional position (x,y,z) 'into' the unit cell,
        i.e.: (x,y,z)->(x',y',z') such that x,y,z in [0,1( """
        pos = np.array(fracpos)
        assert (len(pos) == 3)
        for i in range(3):
            if pos[i]<0:
                while pos[i]<0:
                    pos[i] += 1
            if pos[i]>=1:
                while pos[i]>=1:
                    pos[i] -= 1
        return pos

    def bringFractionalPositionIntoCenteredCell(self, fracpos):
        """Brings a fractional position (x,y,z) 'into' the unit cell,
        i.e.: (x,y,z)->(x',y',z') such that x,y,z in [0,1( """
        pos = np.array(fracpos)
        assert (len(pos) == 3)
        for i in range(3):
            if pos[i]<-0.5:
                while pos[i]<-0.5:
                    pos[i] += 1
            if pos[i]>=0.5:
                while pos[i]>=0.5:
                    pos[i] -= 1
        return pos



    def cartesianPositionInLattice(self, siteId, latticeVector):
        """Returns the cartesian position vector from the origin
        ( fractional coordinates [0,0,0] in unit cell [0,0,0]),
        for a Site corresponding to 'siteId',
        in the unit cell corresponding to latticeVector
        (triplets of coordinates in terms of cellvectors), 
        defining which unit cell in the lattice.
        """
        try:
            posincell = self.getCartesianPosition(siteId)
        except KeyError: raise KeyError( 'Invalid site Id')
        pos = np.array(posincell) + np.dot(latticeVector, self._cellvectors)
        return pos


    def getMonkhorstPackGrid(self, size, shift=(0,0,0)):
        """Returns a Monkhorst-Pack grid of order size[0]*size[1]*size[2],
        scaled by the reciprocal space unit cell.
        The shift is an optional vector shift to all points in the grid."""

        recipvectors = 2 * np.pi * la.inv(np.transpose(self._cellvectors))
        frackpts = MonkhorstPack(size)
        frackpts += np.array(shift)
        # this applies scaling of MP grid by reciprocal cell vectors:
        # (equivalent of frac*vectors[0]+frac*vectors[1]+frac*vectors[2]
        kpts = frackpts*recipvectors.sum(0)
        kpts.shape=(size[0], size[1], size[2], 3)
        return kpts

    def getFracMonkhorstPackGrid(self, size, shift=(0,0,0)):
        """Returns a Monkhorst-Pack grid of order size[0]*size[1]*size[2],
        in fractional coordinates of the reciprocal space unit cell.
        The shift is an optional vector shift to all points in the grid."""

        recipvectors = 2 * np.pi * la.inv(np.transpose(self._cellvectors))
        frackpts = MonkhorstPack(size)
        frackpts += np.array(shift)
        frackpts.shape=(size[0], size[1], size[2], 3)
        return frackpts
        

    pass # end of UnitCell

##########################################################


def create_unitcell( cellvectors, atomList, positionList):
    """Helper function to create a unit cell."""
    rt = UnitCell( cellvectors = cellvectors )
    for a, p in zip( atomList, positionList ):
        site = Site(p,a)
        rt.addSite(site, '')
    return rt


# Here are some tests:

def uc_test1():
    print ("\n*** test1 ***")
    uc = UnitCell( )
    at1=Atom(symbol='Fe', mass=57) ; pos1=(0.0,0.0,0.0)
    at2=Atom(symbol='Al') ; pos2=(0.5,0.5,0.5)
    
    site1 = Site(pos1, at1)
    site2 = Site(pos2, at2)

    uc.addAtom( at1, pos1, "Fe1" )
    uc.addAtom( at2, pos2, "Al1" )
    for site in uc:
       print ("\n position %s \n %s" % (site.getPosition(), site.getAtom()))
       continue
    return


def uc_test2():
    print ("\n*** test2 ***")
    cellvectors = [ (1,0,0), (0,1,0), (0,0,1) ]
    uc = create_unitcell( cellvectors, [Atom(symbol='Fe'), Atom(symbol='Al')], [ (0,0,0), (0.5,0.5,0.5) ] )
    print (uc)
    return

def uc_test3():
    uc = UnitCell( )

    at1 = Atom(symbol='Fe', mass=57)
    at2 = Atom(symbol='Al')
    at3 = Atom(symbol="Zr")

    site1 = Site((0,0,0), at1)
    site2 = Site((0.5,0.5,0.5), at2)
    site3 = Site((0.5, 0.5, 0.0), at3)
    site4 = Site((0.5, 0.0, 0.5), at3)
    site5 = Site((0.0, 0.5, 0.5), at3)
    
    uc.addSite(site1, "Fe1" )
    uc.addSite(site2, "Al1" )
    uc.addSite(site3, "")
    uc.addSite(site4, "")
    uc.addSite(site5, "")

    print ("\n Original unit cell with 3 equivalent Zr atoms:\n")

    for key in uc._siteIds.keys():
        print (key, uc._siteIds[key])

    uc.getSiteFromId("Id0").getAtom().magneticMoment=1.2

    print ("\n Modified unit cell with changed Zr magnetic moment:\n")

    for key in uc._siteIds.keys():
        print (key, uc._siteIds[key])
    return
    
def uc_test4():
    cellvectors = [ (1,0,0), (0,1,0), (0,0,1) ]
    uc = UnitCell(cellvectors=cellvectors)
    assert (uc.getCellVectors() == cellvectors).all()
    return

def test():
    uc_test1()
    uc_test2()
    #uc_test3()
    uc_test4()
    return


if __name__ == "__main__": test()

