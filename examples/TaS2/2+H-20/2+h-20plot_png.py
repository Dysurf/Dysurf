import  numpy as N
import pylab as pl
import matplotlib.pyplot as plt
import math
import sys
plt.style.use('D:/E/E/BIT_conference/experiment/CsBiBr/CsBiBr_kc-small')

fname="SQE_300K"

data=N.loadtxt(fname+".dat")
data = 15.957691216057308*data
fig=plt.figure()
ax=fig.gca()
plt.ylabel('Energy (meV)')
plt.text(0,15.2,"$\mathregular{\Gamma}$")
plt.text(0.4,15.2," M ")
plt.text(0.9,15.2,"$\mathregular{\Gamma}$")
plt.text(1.4,15.2," M ")
plt.text(1.9,15.2,"$\mathregular{\Gamma}$")
plt.xlabel("2+H -2 0 (r.l.u.)")

cax=plt.imshow(pl.rot90(pl.log(data)), interpolation='bicubic', extent=[0.0,2.0,0,15], aspect = 'auto')

plt.set_cmap('jet')

a=7
plt.clim(a,a+6)      ###a: 2..4
ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.yaxis.set_major_locator(plt.MultipleLocator(2))

# plt.colorbar()
plt.savefig(fname+".png",bbox_inches='tight',dpi=600)
