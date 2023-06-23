import numpy as np
import pylab as pl
import os

# Qbh, sma = 1e6, 100.0   # fixed parameters
Qbh, sma = 1e6, 70.0   # fixed parameters

# these need to be the same as "run_Kentr_rhoc.py"
# Kentr_list = [0.2, 0.3, 0.4, 0.5]
# rhocmin, rhocmax = 10, 50
# Nrhoc = 24

Kentr_list = [0.193]
rhocmin, rhocmax = 45, 60
Nrhoc = 24

dir_main = '/Users/wenbinlu/PycharmProjects/Roche_tidal_equilibrium/data_figs/'
color_list = ['tomato', 'royalblue', 'olive', 'orange']
labels = ['qstar', 'xL1', 'PhitotL1', 'xsurf', 'Phitotsurf']
kplt = 4   # which property to be plotted

rhocarr = np.linspace(rhocmin, rhocmax, Nrhoc, endpoint=True)
NK = len(Kentr_list)
status = np.zeros((NK, Nrhoc), dtype=bool)
prop = np.zeros((NK, Nrhoc, 5), dtype=float)   # five properties of the equilibrium solution
for i in range(NK):
    Kentr = Kentr_list[i]
    for j in range(Nrhoc):
        rhoc = rhocarr[j]
        savedir = dir_main + 'Kentr%.3f/rhoc%.3f/' % (Kentr, rhoc)
        fname = savedir + 'output.txt'
        with open(fname, 'r') as f:
            row_new = f.readline()
            row1, row2 = '', ''
            while len(row_new) > 0:
                if 'equilibrium result' in row_new:
                    row1 = row_new
                if 'FINAL' in row_new:
                    row2 = row_new
                row_new = f.readline()
        row1 = row1.replace('equilibrium result: ', '')
        for lab in labels:
            row1 = row1.replace(lab+'=', '').replace(' ', '').strip()
        # print(Kentr, rhoc)
        # print(row1.split(','))
        prop[i, j, :] = [float(item) for item in row1.split(',')]
        prop[i, j, 2] += 1.5*Qbh/sma   # remove the zeroth-order potential term
        prop[i, j, 4] += 1.5*Qbh/sma
        # print(prop[i, j, :])
        if 'detached' in row2:
            status[i, j] = False
        else:
            status[i, j] = True
        # print(status[i, j])

fig, ax = pl.subplots(1, 1, figsize=(13, 9), sharex='all')
for i in range(NK):
    for j in range(Nrhoc):
        if status[i, j]:
            marker = 'x'
        else:
            marker = 'o'
        leg_label = ''
        if j == Nrhoc-1:
            leg_label = 'K=%g' % Kentr_list[i]
        ax.scatter(rhocarr[j], prop[i, j, kplt], s=40, marker=marker,
                   color=color_list[i], label=leg_label)

# if kplt == 3:
#     ax.set_xscale('log')
#     ax.set_yscale('log')

ax.set_xlabel('rhoc')
ax.set_ylabel(labels[kplt])
ax.legend(loc=2, prop={'size': 25}, fancybox=True, framealpha=0.5)
pl.subplots_adjust(bottom=0.115, left=0.12, top=0.98, right=0.98)
fig.savefig(dir_main + labels[kplt] + '.png', dpi=300)

