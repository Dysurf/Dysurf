#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 12})
#matplotlib.rcParams.update({'font.family': "Times"})
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

sqefile = "h04.dat"
bandfile = "omega.dat"
sqe = np.loadtxt(sqefile)
bands = np.loadtxt(bandfile)
qpoints = bands[:,0]
omegas = bands[:,1:]

slice = np.loadtxt("h04.dat")
#print(slice1.shape)
#print(slice2.shape)

ne = 100
dE = 0.5
array_x = qpoints
array_y = np.arange(0, (ne+1)*dE, dE)
color = np.transpose(np.log(sqe))
interpo = interp2d(array_x, array_y, color, kind='cubic')
array_fx = np.linspace(0, qpoints[-1], 400)
array_fy = np.linspace(0, ne*dE, 600)
fcolor = interpo(array_fx, array_fy)
Xgrid, Ygrid = np.meshgrid(array_fx, array_fy, sparse=False, indexing='xy')

fig = plt.figure()
gs = gridspec.GridSpec(5, 2)       # width_ratios=[4,1.5],wspace=.40,hspace=.35

slice_x = np.arange(-12.5, 38, 0.5)
ax = plt.subplot(gs[:2, :])
ax.plot(slice_x, slice[10], color='k', linewidth=2.0, linestyle='-', zorder=1, label=r"$\mathbf{q}$ = 0.05")
ax.plot(slice_x, slice[20], color='r', linewidth=2.0, linestyle='-', zorder=1, label=r"$\mathbf{q}$ = 0.1")
ax.plot(slice_x, slice[40], color='g', linewidth=2.0, linestyle='-', zorder=1, label=r"$\mathbf{q}$ = 0.2")
ax.plot(slice_x, slice[60], color='b', linewidth=2.0, linestyle='-', zorder=1, label=r"$\mathbf{q}$ = 0.3")
ax.plot(slice_x, slice[80], color='m', linewidth=2.0, linestyle='-', zorder=1, label=r"$\mathbf{q}$ = 0.4")
ax.plot(slice_x, slice[99], color='brown', linewidth=2.0, linestyle='-', zorder=1, label=r"$\mathbf{q}$ = 0.5")
ax.set_xlim([-12.5, 38])
#ax.set_xticks(np.linspace(-12.5, 38, 5), minor=False)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(axis='x', which='both', color='k', direction='in')
ax.set_xlabel('Energy (meV)', fontsize=12)
#ax.set_ylim([0, 500])
#ax.set_yticks(np.linspace(0, 500, 6), minor=False)
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(axis='y', which='both', color='k', direction='in')
ax.set_ylabel('Intensity (a.u.)', fontsize=12)
ax.legend(fontsize=12,loc=1,frameon=False)
#ax.text(-1.5,507,r'$\mathbf{a)}$', fontsize=14)

fig.set_size_inches(6*1,6*2)
plt.savefig('slice.eps', bbox_inches="tight", transparent=True)   # dpi=300
#plt.show()


