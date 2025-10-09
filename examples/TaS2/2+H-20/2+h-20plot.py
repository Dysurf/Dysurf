import  numpy as N
import pylab as pl
import matplotlib.pyplot as plt
import math
import sys

fname="SQE_300K"

data=N.loadtxt(fname+".dat")
data = 15.957691216057308*data
plt.interactive(True)

plt.ylabel('Energy (meV)')
plt.text(0,15.2," $\Gamma$ ")
plt.text(0.4,15.2," M ")
plt.text(0.9,15.2," $\Gamma$")
plt.text(1.4,15.2," M ")
plt.text(1.9,15.2," $\Gamma$ ")
plt.text(0.7,-2,"[2+H -2 0](r.l.u.)")

cax=plt.imshow(pl.rot90(pl.log(data)), interpolation='bicubic', extent=[0.0,2.0,0,15], aspect = 'auto')
#                                                                                 #width: aspect
#plt.title("Ta (T=300K)")
plt.set_cmap('jet')

a=7
plt.clim(a,a+6)      ###a: 2..4

plt.colorbar()
# plt.show()
plt.savefig(fname+".eps",bbox_inches='tight')
