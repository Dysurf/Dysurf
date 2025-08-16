# uncompyle6 version 3.7.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.8.2 | packaged by conda-forge | (default, Apr 24 2020, 08:20:52) 
# [GCC 7.3.0]
# Embedded file name: /home/ho8/work/Research/SnSe/SQE/new-CNCS/LDA/RMSD.py
# Compiled at: 2014-10-04 02:08:13


class ReadRootMeanSquareDisplacements:

    def __init__(self, natom):
        self._name = ' Read Root MeanSquare Displacement from a txt file'
        self.displacement = None
        self.natom = natom
        return

    def readrmsd(self, filename):
        infile = open(filename, 'r')
        lines = infile.readlines()
        infile.close()
        first_line = lines[0]
        displacement = []
        i = 0
        st = []
        displacements = []
        for line in lines[1:self.natom + 1]:
            i = i + 1
            st = lines[i].split()
            displacement = [ float(s) for s in st ]
            displacements.append(displacement)

        self.displacements = displacements
# okay decompiling RMSD.pyc
