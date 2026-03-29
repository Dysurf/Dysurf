import  numpy as N
import pylab as pl
import matplotlib.pyplot as pl1
import math
import sys
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Times New Roman'

fname="SQE_14K"

data=N.loadtxt(fname+".dat")
data = 15.957691216057308*data
pl1.interactive(True)

pl1.figure(figsize=(4,3))
pl1.text(-3.1,51,"$\mathregular{\Gamma}$",fontsize=13)
pl1.text(-2.6,51,"M",fontsize=13)
pl1.text(-2.1,51,"$\mathregular{\Gamma}$",fontsize=13)
pl1.text(-1.6,51,"M",fontsize=13)
pl1.text(-1.1,51,"$\mathregular{\Gamma}$",fontsize=13)
pl1.text(-2.5,-7.45,"[H-H2] (r.l.u.)",fontsize=13)
pl1.text(-2.95,46.5,"14K",fontsize=13,color='white')

cax=pl1.imshow(pl.rot90(pl.log(data)), interpolation='bicubic', extent=[-3.0,-1.0,0,50], aspect = 0.05)
pl1.set_cmap('jet')

pl1.clim(3,9)
pl1.yticks([])

pl1.colorbar()
# pl1.show()
pl1.savefig(fname+".png",bbox_inches='tight')

