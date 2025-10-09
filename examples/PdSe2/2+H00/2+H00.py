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
plt.text(0,16.2,"$\mathregular{\Gamma}$")
plt.text(0.48,16.2," X ")
plt.xlabel("2+H 0 0 (r.l.u.)")

cax=plt.imshow(pl.rot90(pl.log(data)), interpolation='bicubic', extent=[0.0,0.5,0,16], aspect = 'auto')

plt.set_cmap('jet')

a=7
plt.clim(a,a+6)      ###a: 7..13
ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
ax.yaxis.set_major_locator(plt.MultipleLocator(2))

# plt.colorbar()
plt.savefig(fname+".eps",bbox_inches='tight',dpi=600)
