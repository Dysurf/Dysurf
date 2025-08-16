# [001] and [110] plane
import sys
import numpy as N
import math
import cmath
from scipy import signal
import pylab as pl
import matplotlib.pyplot as pl1
from crystal.UnitCell import *
from crystal.SimpleAtom import *
import kernelGenerator2
from idf.FractionalQs import read as readQs
import phonopy.file_IO as file_IO
from phonopy.structure.atoms import Atoms as PhonopyAtoms
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from RMSD import *
from ResFuncsCNCS import *
from kernelGenerator2.sqeForExistedPhononsForGivenQ import SQECalculatorForExistingPhononsForGivenQ
from datetime import datetime
import rotation
import pickle

####add
#from phonopy.structure.symmetry import symmetrize_borns_and_epsilon
from phonopy.structure.symmetry import Symmetry
from phonopy.units import Hartree, Bohr
####
starttime=datetime.now()

#==============================================
# Read the root mean square displacment (RMSD)
#==============================================
" eneter the number of atoms in the primitive unit cell"
natom=4
RMSDFile=ReadRootMeanSquareDisplacements(natom)
RMSDFile.readrmsd('RMSD.txt')
displacements=RMSDFile.displacements
print(displacements)
#==================================================


force_constants=file_IO.parse_FORCE_CONSTANTS('FORCE_CONSTANTS')
bulk = read_vasp("POSCAR")
lattice_vector=bulk.cell
print("lattice", lattice_vector)
phonon=Phonopy(bulk,[[4,0,0], [0,4,0], [0,0,3]])

####################################################################
# For different sysmtes, you need to change something here !!
#####################################################################

primitive_cell=N.array([[3.1900000572,0,0],[-1.5950000286,2.7626210876,0],[0,0,5.1900000572]])
inv_primitive_cell=N.linalg.inv(primitive_cell)*2*N.pi


#print "delta=",delta

ne = 500      # Number of energy  Emax=ne*deltae
deltae = 0.1  # delta of energy


uc = UnitCell( )
uc.setCellVectors([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
#uc.setCellVectors(primitive_cell.transpose())
at1=Atom(symbol='Ga',Z=7.288,  mass=69.72)
at1.coh_sct_length=7.288
at1.incoh_sct_length=0.0
#site1 = Site(pos1, at1)


at2=Atom(symbol='N',Z=9.36,  mass=14.007)
at2.coh_sct_length=9.36
at2.incoh_sct_length=0.0
#site2 = Site(pos3, at3)


pos1=( 0.333333357,  0.666666714,  0.000000000 )
pos2=( 0.666666620,  0.333333314,  0.500000000 )
pos3=( 0.333333357,  0.666666714,  0.376500004 )
pos4=( 0.666666620,  0.333333314,  0.876500004 )


uc.addAtom( at1, pos1, "Ga" )
uc.addAtom( at1, pos2, "Ga" )
uc.addAtom( at2, pos3, "N" )
uc.addAtom( at2, pos4, "N" )


#uc.addAtom( at1, pos2, "Sn" )
#uc.addAtom( at1, pos3, "Sn" )
#uc.addAtom( at1, pos4, "Sn" )
#uc.addAtom( at2, pos5, "Se" )
#uc.addAtom( at2, pos6, "Se" )
#uc.addAtom( at2, pos7, "Se" )
#uc.addAtom( at2, pos8, "Se" )

##  You don't need to change anything below here for the calculation.
##  But you may need to change something in Plot section to make it look better
##  Also, if you need to "shift" the bin box, you also need to change some thing
##  in the section of "Bin box calculaiton"
##  In principle, you should shift bin box easily. But I haven't thought about it yet.
###################################################################################

symmetry = phonon.get_symmetry()
##inceased by myself
phonon.generate_displacements(distance=0.01)
##
supercells = phonon.get_supercells_with_displacements()
born=[[[5.9746620,0.0000000, 0.0000000],
[0.0000000,5.3734410,0.0000000],
[0.0000000,0.0000000,4.9273880]],
[[ 1.1377200,0.0000000,0.0000300],
[0.0000000,1.1703100 , 0.0000400],
[0.0000000, 0.0000000,1.1487100]], 
[[-7.1987600 , 0.0000000 , 0],
[-0.0000000,-1.3200500,0.0000900],  
[-0.0000000,-0.0000000 ,-1.5118200]],
[[-6.9099500,-0.0000000, 0],
[0.0000000,-1.4667500 ,0],
[-0.0000000, 0.0000000 ,-1.3998400]],
[[-1.5806900, 0, 0],
[0.0000000 ,-3.7893500, 2.2268600],
[0.0000000, 2.0581500, -3.2718000]],
[[-1.4755300,0, -0],
[0.0000000, -3.9025900, 2.3767400], 
[ -0.0000000 , 2.1369600,-3.2747800]],   
[[1.1109400, 0.0000000, 0.00 ],       
[0.0000000 , 1.2196300 , 0.0000],       
[ -0.0000000 , 0.0000000 ,1.1317500]],    
[[9.2142000, -0.0000000 , 0.000 ],    
[0.0000000, 7.8616300, 0.0000],   
[-0.0000000 ,-0.0000000, 6.9705400]],    
[[8.7651900,-0.0000000 ,0.00000   ], 
[ 0.0000000, 7.9302600, 0.0000],     
[0.0000000 , 0.0000000 , 6.7638700]]]
 
 
epsilon=[[6.511285, 0.000000, 0.000000],
[0.000000, 6.511285, 0.000000],
[0.000000, 0.000000, 6.511285]]

factors=14.400
#
phonon.set_force_constants(force_constants)
#
#phonon.set_nac_params({'born':born,'factor':factors,'dielectric':epsilon})
#phonon.produce_force_constants()
#phonon.set_post_process(primitive_cell,force_constants=force_constants,is_nac=True)

print("Space group:", symmetry.get_international_table())

VasptoTHz=15.633302
THz2meV=4.1357
bands=[]

# Generate q points in bin box (ql,qt1,qt2) and convert it to (qx,qy,qz) for the
# band struture calculation
fname="SQE_GaN"
nqlong=400 
nqperpen=10 # 0 --6
#
delta_h=0.005
delta_k=0.01
delta_l=0.01

nqh= nqlong
nqk= nqperpen
nql= nqperpen


####################
x1=-3
y1=3
z1=2

nqtot=0              # total number of q points

for ih in range(nqh):
#   qh=delta_h*(ih-nqh/2)
   qh=delta_h*ih
   for ik in range(nqk):
       qk=delta_k*(ik-nqk/2)
       for il in range(nql):
           ql=delta_l*(il-nql/2)
#########
           qx=x1+qh+qk
           qy=y1-qh+qk
           qz=z1+ql
#           R=rotation.rotation([h,k,l],lattice_vector)  # Rotation matrix
#           q=N.dot(R,[ql,qt1,qt2])+[h1,k1,l1]           # transfer coordinate in bin box (ql,qt1,qt2) to (qx,qy,qz)
#           q=N.array([N.dot(primitive_cell,q)])
           q=N.array([[qx,qy,qz]])
           nqtot=nqtot+1
           bands.append(q)
print("Total of qpoint", nqtot)

bands=N.array(bands)

phonon.set_band_structure(bands,is_eigenvectors=True)
output=phonon.get_band_structure()
distances = output[1]
eigval = output[2]
eigvecs = output[3]
qpoints = output[0]

energies = N.vstack(eigval)*THz2meV

energies[N.where(energies < 0)] = 0
energies[N.where(energies < 1e-3)] += 1e-3

qvecs = N.vstack(qpoints)
all_dists = N.hstack(distances)
eigvecs1 = N.concatenate(eigvecs)

pol=[x.T for x in eigvecs1]
pol=N.array(pol)
pol=pol.reshape(nqtot,3*natom,natom,3)

#convert to Cartesian coordinates
for i in range(nqtot):
    qvecs[i]=N.dot(inv_primitive_cell,qvecs[i])

#===================================

sqecalc = kernelGenerator2.SqeCalculator.SqeCalculator(uc, kpoints=qvecs,energies=energies,polarizations=pol,displacements=displacements)
sqeExistedPhononsForGivenQ=SQECalculatorForExistingPhononsForGivenQ(sqecalc)
qt=qvecs
g = N.zeros(ne)
evalues=N.arange(0,(ne-1)*deltae+0.01,deltae)
ecenter=N.zeros(ne)
gaussRes=GaussianResolution()

#########################  Bin box calculation ######################################## 
BinnedSQE=N.zeros((nqk,nql,nqh,ne))
#BinnedSQE=N.zeros((nqh,nql,ne))

ntotal=0              # total number of q points

for ih in range(nqh):
   for ik in range(nqk):
       for il in range(nql):
            output=sqeExistedPhononsForGivenQ.calcSqeCohforExistedPhononsforGivenQ(qt[ntotal])
            for j in range(len(output)) :
                EIndex=int(round(output[j][0]/deltae))
                Convolution = Gauss12(evalues,qt[ntotal],output[j][0],output[j][1],gaussRes.sigmae12,0.01)  # 0.025
                BinnedSQE[ik,il,ih,:]+=output[j][2]*Convolution

            ntotal=ntotal+1
            if (ntotal%1000==0) :#print status only every 1000 points
               print(ntotal,' of ',nqtot,' ',ntotal*100.0/nqtot,' percent')

BinnedSQE[N.where(BinnedSQE<1e-3)] += 1e-3

######################### End of bin box calculation ########################################

#sum2=BinnedSQE[nt/2,nt/2,:,:]

sum3=sum(sum(BinnedSQE[:,:,:,:],0),0)

#N.savetxt(fname+"_half.dat",sum2)

#N.savetxt(fname+"_half.dat",sum2)
N.savetxt(fname+'.dat',sum3)

#N.save(fname,BinnedSQE)
#N.savetxt(fname+'.dat',BinnedSQE)

#pickle.dump(BinnedSQE,"test")
stoptime=datetime.now()
deltatime=stoptime-starttime

print(deltatime)

#import matplotlib.pyplot as pl1
#fname="H02-202-00L-LOTO"
#



#clim1=sys.argv[1]
#clim2=sys.argv[2]

#data=N.load(fname+".npy")
#data=sum(sum(data[:,:,:,:],0),0)
data=N.loadtxt(fname+".dat")

pl1.interactive(True)
pl1.text(0,-3.5,"(2 0 2)-->>(2 0 2+k)",fontsize=12)
pl1.text(0,-5.5,"Gamma-->Z",fontsize=12)
pl1.ylabel('Energy (meV)',fontsize=17)
cax=pl1.imshow(pl.rot90(pl.log(data)), interpolation='bicubic', extent=[-3,-1,0,50], aspect = 0.02)
#                                                                                 #width: aspect
#pl1.title("Ta (T=300K)")
pl1.set_cmap('jet')

pl1.title(fname+"(T=300K)",fontsize=17)
 #pl1.clim(5,11)
#pl1.clim(7,12)
a=10
pl1.clim(a,a+4)      ###a: 2..4

pl1.colorbar()
pl1.show()
#pl1.savefig(fname+".eps")
#raw_input("input to end this one")

