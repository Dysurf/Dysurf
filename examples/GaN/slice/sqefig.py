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

sqefile1 = "SQE_300K_LA.dat"
bandfile1 = "omega_LA.dat"
sqefile2 = "SQE_300K_TA.dat"
bandfile2 = "omega_TA.dat"
sqe1 = np.loadtxt(sqefile1)
bands1 = np.loadtxt(bandfile1)
qpoints1 = bands1[:,0]
omegas1 = bands1[:,1:]
sqe2 = np.loadtxt(sqefile2)
bands2 = np.loadtxt(bandfile2)
qpoints2 = bands2[:,0]
omegas2 = bands2[:,1:]

slice1 = np.loadtxt("SQE_300K1.dat")
slice2 = np.loadtxt("SQE_300K2.dat")
#print(slice1.shape)
#print(slice2.shape)

ne = 300
dE = 0.05
array_x1 = qpoints1
array_x2 = qpoints2
array_y = np.arange(0, (ne+1)*dE, dE)
color1 = np.transpose(np.log(sqe1))
interpo1 = interp2d(array_x1, array_y, color1, kind='cubic')
color2 = np.transpose(np.log(sqe2))
interpo2 = interp2d(array_x1, array_y, color2, kind='cubic')
array_fx1 = np.linspace(0, qpoints1[-1], 400)
array_fx2 = np.linspace(0, qpoints2[-1], 400)
array_fy = np.linspace(0, ne*dE, 600)
fcolor1 = interpo1(array_fx1, array_fy)
fcolor2 = interpo2(array_fx1, array_fy)
Xgrid1, Ygrid = np.meshgrid(array_fx1, array_fy, sparse=False, indexing='xy')
Xgrid2, Ygrid = np.meshgrid(array_fx2, array_fy, sparse=False, indexing='xy')

fig = plt.figure()
gs = gridspec.GridSpec(5, 2)       # width_ratios=[4,1.5],wspace=.40,hspace=.35

slice_x = np.arange(0, (300+1)*0.05, 0.05)
ax = plt.subplot(gs[:2, :])
ax.plot(slice_x, slice1[33]+5, color='k', linewidth=2.0, linestyle='-', zorder=1, label=r"$\mathbf{Q}$ = (1.1,1.1,0)")
ax.plot(slice_x, slice2[33]+20, color='r', linewidth=2.0, linestyle='-', zorder=2, label=r"$\mathbf{Q}$ = (1.2,0.9,0)")
ax.set_xlim([0.0, 12])
ax.set_xticks(np.linspace(0, 12, 7), minor=False)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(axis='x', which='both', color='k', direction='in')
ax.set_xlabel('Energy (meV)', fontsize=12)
ax.set_ylim([0, 200])
ax.set_yticks(np.linspace(0, 200, 6), minor=False)
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(axis='y', which='both', color='k', direction='in')
ax.set_ylabel('Intensity (a.u.)', fontsize=12)
ax.legend(fontsize=12,loc=1,frameon=False)
ax.text(-1.5,207,r'$\mathbf{a)}$', fontsize=14)

ax = plt.subplot(gs[2:, 0])
ax.pcolormesh(Xgrid1, Ygrid, fcolor1, cmap='jet', vmin=0, vmax=5, zorder=0, shading='gouraud')
ax.plot(qpoints1, omegas1, color='w', linewidth=1.0, linestyle='-', zorder=1)
ax.set_xlim([0.0, qpoints1[-1]])
ax.set_xticks([0.0, qpoints1[-1]/2, qpoints1[-1]], minor=False)
ax.set_xticklabels([r'$\Gamma$',r'$[1+\xi,1+\xi,0]$','K'])
ax.tick_params(axis='x', which='both', color='w', direction='in')
#ax.set_xlabel('Wavevector', fontsize=14)
ax.set_ylim([0.0, 12.5])
ax.set_yticks(np.linspace(0, 12, 5), minor=False)
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(axis='y', which='both', color='w', direction='in')
ax.set_ylabel('Energy (meV)', fontsize=12)
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.text(-0.30,12.0,r'$\mathbf{b)}$', fontsize=14)

ax = plt.subplot(gs[2:, 1])
ax.pcolormesh(Xgrid2, Ygrid, fcolor2, cmap='jet', vmin=0, vmax=5, zorder=0, shading='gouraud')
ax.plot(qpoints2, omegas2, color='w', linewidth=1.0, linestyle='-', zorder=1)
ax.set_xlim([0.0, qpoints2[-1]])
ax.set_xticks([0.0, qpoints2[-1]/2, qpoints2[-1]], minor=False)
ax.set_xticklabels([r'$\Gamma$',r'$[\xi,\xi,4]$','K'])
ax.tick_params(axis='x', which='both', color='w', direction='in')
#ax.set_xlabel('Wavevector', fontsize=14)
ax.set_ylim([0.0, 12.5])
ax.set_yticks(np.linspace(0, 12, 5), minor=False)
ax.set_yticklabels([])
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(axis='y', which='both', color='w', direction='in')
#ax.set_ylabel('Energy (meV)', fontsize=12)
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.text(-0.15,12.0,r'$\mathbf{c)}$', fontsize=14)
#plt.tight_layout()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
fig.set_size_inches(6*1.2, 6*1.6)
plt.savefig('Fig.pdf', bbox_inches="tight", transparent=True)   # dpi=300
#plt.show()


