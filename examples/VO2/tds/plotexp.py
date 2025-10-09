import  numpy as N
import pylab as pl
import matplotlib.pyplot as pl1
import math
import sys


x = N.genfromtxt('M.txt', dtype=float, delimiter=' ', usecols = (0))
y = N.genfromtxt('M.txt', dtype=float, delimiter=' ', usecols = (1))
z = N.genfromtxt('M.txt', dtype=float, delimiter=' ', usecols = (2))
z[N.where(z<0)]=0
print ("zmax=",max(z), "zmin=",min(z))
z=z/max(z)
print ("zmax=",max(z), "zmin=",min(z))

Xnew=N.reshape(x,(142,154))   # "1.txt"
Ynew=N.reshape(y,(142,154))
Znew=N.reshape(z,(142,154))

#clim1=sys.argv[1]
#clim2=sys.argv[2]



#x = N.genfromtxt('2.txt', dtype=float, delimiter=' ', usecols = (0))
#y = N.genfromtxt('2.txt', dtype=float, delimiter=' ', usecols = (1))
#z = N.genfromtxt('2.txt', dtype=float, delimiter=' ', usecols = (2))
#Xnew=N.reshape(x,(112,157))   # "1.txt"
#Ynew=N.reshape(y,(112,157))
#Znew=N.reshape(z,(112,157))



#pl.pcolormesh(Xnew,Ynew,Znew)
#pl.show()



pl1.interactive(True)
cax=pl1.imshow(pl.rot90(Znew), interpolation='bicubic' ,aspect=0.68) #, extent=[50,70,0,90])# , aspect = ratio)
#cax=pl1.imshow(Znew, interpolation='bicubic' ,aspect=0.625) #, extent=[0,1,0,1])# , aspect = ratio)

pl1.xlim(47,86)
pl1.ylim(125,30)

color=sys.argv[1]
pl1.set_cmap(color)

#pl1.title("VO2 (T=425 K), R phase")
#pl1.title("VO2"+"("+sys.argv[2]+","+sys.argv[3]+")")

#pl1.clim(clim1,clim2)

pl1.clim(0., 0.035) # better


#climlim1=sys.argv[1]
#im1=sys.argv[1]
#pl1.clim(0,0.01)   # works well for production

vmin, vmax = pl1.gci().get_clim()

print ("vmin=",vmin,"vmax=",vmax)


#pl1.xlabel("[1 -1 0]")
#pl1.ylabel("[0 0 1]")


#pl1.text(50, 85, "(-5 4 3)",color="white")
#pl1.text(90, 85, "(-4 5 3)",color="white")
#pl1.text(70, 115, "(-5 5 2)",color="black")

#pl1.xticks([],(""))
#pl1.yticks([],(""))


#pl1.xticks([40,68],("-5","-4"))
#pl1.yticks([85,50,10],("2","3","4")) #,rotation=90)

#pl1.xticks([16,56,96,136],("-6","-5","-4","-3"))
#pl1.yticks([8,48,88],("4","3","2")) #,rotation=90)
#pl1.text(10, 30, "(-5 4 3)",rotation=90,color="white")
#pl1.text(4, 40, "(-5 4 3)")
#pl1.text(56, 40, "(-4 5 3)")
#pl1.text(30, 16, "(-4 4 4)")
#pl1.text(30, 66, "(-5 5 2)")

#pl1.subplots_adjust(left=0.,right=0.,bottom=0.,top=0.,wspace=0.,hspace=0.)



pl1.xticks([],(""))
pl1.yticks([],(""))


#pl1.colorbar()

#pl1.colorbar()
#pl1.show()
pl1.savefig("exp.eps",bbox_inches='tight')
#raw_input("input to end this one")



