#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 12})
#matplotlib.rcParams.update({'font.family': "Times"})
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np

matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

sqefile = "SQE_300K.dat"
bandfile = "omega.dat"
sqe = np.loadtxt(sqefile)
bands = np.loadtxt(bandfile)
qpoints = bands[:,0]
omegas = bands[:,1:]

ne = 300
dE = 0.05
array_x = qpoints
array_y = np.arange(0, (ne+1)*dE, dE)
color = np.transpose(np.log(sqe))
interpo = interp2d(array_x, array_y, color, kind='cubic')
array_fx = np.linspace(0, qpoints[-1], 1000)
array_fy = np.linspace(0, ne*dE, 1000)
fcolor = interpo(array_fx, array_fy)
Xgrid, Ygrid = np.meshgrid(array_fx, array_fy, sparse=False, indexing='xy')

plt.figure()
plt.pcolormesh(Xgrid, Ygrid, fcolor, cmap='jet', zorder=0, shading='gouraud')
plt.plot(qpoints, omegas, color='w', linewidth=1, linestyle='-', zorder=1)
ax = plt.gca()
ax.set_xlim([0.0, qpoints[-1]])
ax.set_xticks([0.0, qpoints[-1]/2, qpoints[-1]], minor=False)
ax.set_xticklabels([r'$\Gamma$','[1+H,1+H,0]','M'])
ax.tick_params(axis='x', which='both', color='w', direction='in')
#ax.set_xlabel('Wavevector', fontsize=14)
ax.set_ylim([0.0, ne*dE])
ax.set_yticks(np.linspace(0, ne*dE, 6), minor=False)
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(axis='y', which='both', color='w', direction='in')
ax.set_ylabel('Energy (meV)', fontsize=12)
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
plt.clim(0, 6)
plt.colorbar()

plt.savefig(sqefile.split('.')[0]+'.pdf', bbox_inches='tight')
#plt.show()


