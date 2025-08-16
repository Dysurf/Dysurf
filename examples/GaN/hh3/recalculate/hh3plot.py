import  numpy as N
import pylab as pl
import matplotlib.pyplot as pl1
import math
import sys

fname="hh3_thick"
#clim1=sys.argv[1]
#clim2=sys.argv[2]

#data=N.load(fname+".npy")
#data=sum(sum(data[:,:,:,:],0),0)
data=N.loadtxt(fname+".dat")
data = 15.957691216057308*data
pl1.interactive(True)
#pl1.xlabel('Q=(0 0 2)-->>(2+k  -2-k  2+k)')
#pl1.text(0,-1.5,"(0 0 0)-->>(h 0 0)",fontsize=12)
#pl1.text(0,-2,"Gamma-->H",fontsize=12)
#pl1.ylabel('Energy (meV)',fontsize=17)
pl1.text(-0.3,51," K ",font='Times New Roman',fontsize=17)
pl1.text(0,51," $\Gamma$ ",font='Times New Roman',fontsize=17)
pl1.text(0.3,51," K ",font='Times New Roman',fontsize=17)
#pl1.text(0.15,17.4,"PdSe2(T=300K)",fontsize=17)
pl1.text(-0.2,-7,"[HH3] (r.l.u.)",fontsize=17)
pl1.text(-0.47,46.5,"300K",fontsize=16,color='white')

cax=pl1.imshow(pl.rot90(pl.log(data)), interpolation='bicubic', extent=[-0.5,0.5,0,50], aspect = 0.025)
#pl1.title("Ta (T=300K)")
pl1.set_cmap('jet')
#pl1.title("CsI(T=0K)",fontsize=12)
#a=3
#pl1.clim(a,a+3)      ###a: 2..4
#pl1.text(1,-2.7,"$H^2$+$K^2$+$L^2$",fontsize=17)

pl1.clim(3,9)
pl1.xticks(fontsize=15)
pl1.yticks([])
#pl1.text(0.0,-2.7,"G",fontsize=17)
#pl1.text(0.5,-2.7,"M",fontsize=17)
#pl1.text(0.97,-2.7,"G",fontsize=17)
#pl1.clim(clim1,clim2)
#pl1.colorbar()
pl1.show()
pl1.savefig(fname+".eps",bbox_inches='tight')
#pl1.savefig(fname+".eps")

#pl1.savefig(fname+".eps")
#pl1.savefig("test.eps")

#raw_input("input to end this one")


