#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

data = np.loadtxt('SQEBIN_300K.dat')
data2=pd.DataFrame(data)
#print(data2)
#data3=data[]
#leng=len(data)
#for i in range(leng):
#  print (data)
data3=data2[data2[2]-0.025<=1e-7]
data3.loc[:,5]=np.log(data3.loc[:,4])
print(data3)
data3.to_csv("log4D_try.csv",index=False,sep=' ',header=0)
###np.savetxt('log4D.dat',data3)

#lll=pd.Series(np.log(data3.loc[:,4].values))
#print(data3.loc[:,4].values())
#print(lll)
#data4=pd.concat([data3,],axis=1)
#print(data4)
#print(pd.Series(np.log(data3.loc[:,4].values)))
#data4.to_csv('log4D.dat',sep='\t')
#  if data[i][3]- 0.5000000000000000E-01<=1e-7:
#    np.savetxt('4D.dat',data[i])
#print("OK")
