#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 13:31:25 2023

@author: egustafs
"""

import pandas 
import numpy as np
import gvar as gv
import lsqfit
import matplotlib.pyplot as plt
import ancillary

import matplotlib as mpl
from  matplotlib import rc
import matplotlib.pyplot as plt


colorsrgb = [(0, 0, 0), (230, 159, 0), (86, 180, 233),
          (0, 158, 115), (240, 228, 66), (0, 114, 178),
          (213, 94, 0), (204, 121, 167)]

fmts = ['o', 'd', 's', '^', '<', '8', '6', 'd', 'x']
colors = [tuple([el / 255 for el in c]) for c in colorsrgb]

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = "Times New Roman"
mpl.rcParams.keys()
# mpl.rcParams["text.fontset"] = "Times New Roman"
mpl.rcParams["mathtext.fontset"] = "stix"


datapoints = gv.load('readoutcorrelations.pkl')[(7, 2)]

for k in range(20):
    plt.figure()
    plt.title(str(k))
    points = datapoints[:, k]
    
    plt.errorbar(np.linspace(1, 30, 30), gv.mean(points), yerr=gv.sdev(points),
                 fmt='.', capsize=5, markersize=10, color=colors[0],
                 label='Observable')
    plt.xlabel('Circuit Number')
    plt.ylabel(r'$\langle O_i(t=16) \rangle$')
    
    wavg = lsqfit.wavg(points)
    
    plt.fill_between([0.5, 30.5], [gv.mean(wavg) + gv.sdev(wavg)] * 2,
                     [gv.mean(wavg) - gv.sdev(wavg)] * 2,
                     facecolor=[colors[1][i] for i in range(3)] + [0.5],
                     hatch='/',
                     label='correlated weighted average')
    
    
    avg = lsqfit.wavg(gv.gvar(gv.mean(points), gv.sdev(points)))
    print(avg)
    plt.fill_between([0.5, 30.5], [gv.mean(avg) + gv.sdev(avg)] * 2,
                     [gv.mean(avg) - gv.sdev(avg)] * 2,
                     facecolor=[colors[2][i] for i in range(3)] + [0.5],
                     hatch='|',
                     label='uncorrelated weighted average')
    
    
    avg = gv.gvar(np.mean(gv.mean(points)), 
                  np.std(gv.mean(points)) / np.sqrt(30))
    print(avg)
    print('='*60)
    plt.fill_between([0.5, 30.5], [gv.mean(avg) + gv.sdev(avg)] * 2,
                     [gv.mean(avg) - gv.sdev(avg)] * 2,
                     facecolor=[colors[3][i] for i in range(3)] + [0.5],
                     hatch='x',
                     label='average')
    
    plt.legend(framealpha=0)
    plt.savefig('correlatedaveragecomparison.pdf')


points = gv.evalcorr(datapoints[:, 2])
x, y = np.meshgrid(np.linspace(1, 30, 30), np.linspace(1, 30, 30))

fig, ax = plt.subplots()
c = ax.pcolormesh(x, y, points, vmin=0, vmax=0.1)
fig.colorbar(c, ax=ax)
ax.set_title('Correlations at equal time for different twirled circuits')
ax.set_xlabel('Twirled Circuit Number')
ax.set_ylabel('Twirled Circuit Number')
fig.savefig('readoutcorrelationsequaltime.pdf')
# lsqfit.wavg()




# points = gv.evalcorr(datapoints[7])

# x, y = np.meshgrid(np.linspace(1, 20, 20), np.linspace(1, 20, 20))

# fig, ax = plt.subplots()
# c = ax.pcolormesh(x, y, points, vmin=0, vmax=0.1)
# fig.colorbar(c, ax=ax)
# ax.set_title('Readout correlations at different time steps')
# ax.set_xlabel('Time step')
# ax.set_ylabel('Time step')
# fig.savefig('readoutcorrelationsacrosstime.pdf')