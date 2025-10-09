import  numpy as N
import pylab as pl
import matplotlib.pyplot as plt
import math
import sys
import matplotlib as mpl
plt.style.use('D:/E/E/BIT_conference/experiment/CsBiBr/CsBiBr_kc-small') # this can be commented out.

fname="SQE_0K"

data=N.loadtxt(fname+".dat")
data = 15.957691216057308*data

fig=plt.figure()
ax=fig.gca()
plt.ylabel('Energy (meV)')
plt.text(0,12.2,"$\mathregular{\Gamma}$")
plt.text(1,12.2,"X")
plt.text(2,12.2,"$\mathregular{\Gamma}$")
plt.text(3,12.2,"X")
plt.text(3.8,12.2,"$\mathregular{\Gamma}$")

cax=plt.imshow(pl.rot90(pl.log(data)), interpolation='bicubic', extent=[0.0,4.0,0,12], aspect = 0.34)

plt.set_cmap('jet')
plt.xlabel("H 1 0 (r.l.u.)")

a=4.5
plt.clim(a,a+3)     
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.yaxis.set_major_locator(plt.MultipleLocator(2))

plt.colorbar()
plt.savefig(fname+".png",bbox_inches='tight',dpi=600)
