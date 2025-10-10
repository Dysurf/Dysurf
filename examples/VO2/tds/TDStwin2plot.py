import  numpy as N
import pylab as pl
import matplotlib.pyplot as pl1
import math
import sys

#print "usage: python plotSQ.py filename clim_min clim_max k_min K_max, eg: SQ1.4 8 12 2.4 3.5"

fname="TDStwin2"

data = N.loadtxt(fname+".dat")
data = 15.957691216057308*data
#data = data*15.3
pl1.interactive(True)
cax=pl1.imshow(pl.rot90(pl.log(data)), interpolation='bicubic' ,aspect=1.5) #, extent=[0,1,0,1])# , aspect = ratio)

pl1.xlim(2,27)
pl1.ylim(27,3)

color= "jet" #sys.argv[1]
pl1.set_cmap(color)

pl1.clim(9,21) # better

pl1.xlabel(r"0 $\mathrm{\xi}$ $\mathrm{\xi}$ direction",fontsize=14)
pl1.ylabel(r"4$\mathrm{\xi}$ $\mathrm{\xi}$ -3$\mathrm{\xi}$ direction",fontsize=14)

li=[3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6, 12.8, 13.0, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.4, 14.6, 14.8, 15.0, 15.2, 15.4, 15.6, 15.8, 16.0, 16.2, 16.4, 16.6, 16.8, 17.0]

liy=[-5.9, -5.8, -5.7, -5.6, -5.5, -5.4, -5.3, -5.2, -5.1, -5.0, -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

pl1.xticks([2,14,27],[li[2],li[14],li[27]],fontsize=14)
pl1.yticks([27,15,3],[liy[3],liy[15],liy[27]],fontsize=14)

pl1.savefig(fname+".eps",bbox_inches='tight')
pl1.savefig(fname+".jpg",dpi=600,bbox_inches='tight')

