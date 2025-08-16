import  numpy as N
import pylab as pl
import matplotlib.pyplot as pl1
import math
import sys

fname="1+h10"
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
pl1.ylabel('Energy (meV)',fontsize=17)
pl1.text(0,15.2," $\Gamma$ ",fontsize=17)
pl1.text(0.4,15.2," M ",fontsize=17)
pl1.text(0.9,15.2," $\Gamma$ ",fontsize=17)
pl1.text(1.4,15.2," M ",fontsize=17)
pl1.text(1.9,15.2," $\Gamma$ ",fontsize=17)
pl1.text(0.7,16.4,"TaS2(T=300K)",fontsize=17)
pl1.text(0.7,-2,"[1+H 1 0](r.l.u.)",fontsize=17)

cax=pl1.imshow(pl.rot90(pl.log(data)), interpolation='bicubic', extent=[0.0,2.0,0,15], aspect = 'auto')
#                                                                                 #width: aspect
#pl1.title("Ta (T=300K)")
pl1.set_cmap('jet')

#pl1.title("CsI(T=0K)",fontsize=12)
#pl1.clim(5,11)
#pl1.clim(7,12)
a=7
pl1.clim(a,a+6)      ###a: 2..4
#pl1.text(1,-2.7,"$H^2$+$K^2$+$L^2$",fontsize=17)

#pl1.clim(8.5,11.5)
#pl1.text(0.0,-2.7,"G",fontsize=17)
#pl1.text(0.5,-2.7,"M",fontsize=17)
#pl1.text(0.97,-2.7,"G",fontsize=17)
#pl1.clim(clim1,clim2)
pl1.colorbar()
pl1.show()
pl1.savefig(fname+".eps",bbox_inches='tight')
#pl1.savefig(fname+".eps")

#pl1.savefig(fname+".eps")
#pl1.savefig("test.eps")

#raw_input("input to end this one")


