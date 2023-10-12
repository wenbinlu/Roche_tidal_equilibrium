import numpy as np
import pylab as pl
import os
from scipy.interpolate import interp1d
from glob import glob

Qbh, sma = 1e6, 150.0   # fixed parameters

# these need to be the same as "run_Kentr_rhoc.py"
# Kentr_list = [0.2, 0.3, 0.4, 0.5]
# rhocmin, rhocmax = 10, 50
# Nrhoc = 24

qstar_target = 0.5

Kentr_list = [0.157]
Nrhoc = 24

dir_main = '/Users/wenbinlu/PycharmProjects/Roche_tidal_equilibrium/data_figs/sma%d/' % sma
color_list = ['tomato', 'royalblue', 'olive', 'orange']
labels = ['qstar', 'xL1', 'PhitotL1', 'xsurf', 'Phitotsurf', 'rhoc']
kplt_list = [4]   # which property to be plotted [0-5]

NK = len(Kentr_list)
status = np.zeros((NK, Nrhoc), dtype=bool)
prop = np.zeros((NK, Nrhoc, 6), dtype=float)   # five properties of the equilibrium solution
for i in range(NK):
    Kentr = Kentr_list[i]
    # print(os.walk(dir_main + 'Kentr%.3f/'))
    # dir_list = [x[0] for x in os.walk(dir_main + 'Kentr%.3f/')]
    dir_list = glob(dir_main + 'Kentr%.3f/*/' % Kentr, recursive=True)
    # print(dir_list)
    for j in range(len(dir_list)):
        savedir = dir_list[j]
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
    rhocarr = prop[i, :, 5]
    xsurfarr = prop[i, :, 3]
    Phitotsurf = prop[i, :, 4]
    qstararr = prop[i, :, 0]
    if max(qstararr) > qstar_target > min(qstararr):
        intp_Phisurf = interp1d(qstararr, Phitotsurf)
        intp_xsurf = interp1d(qstararr, xsurfarr)
        print('Phitotsurf(qstar_target=%.5f)+1.5Qbh/sma=%.5f'
              % (qstar_target, intp_Phisurf(qstar_target)))
        print('xsurf(qstar_target=%.5f)=%.5f'
              % (qstar_target, intp_xsurf(qstar_target)))
    else:
        print('qstar_target outof simulated range: [%.5f, %.5f]' % (min(qstararr), max(qstararr)))
    for j in range(Nrhoc):
        if status[i, j]:
            marker = 'x'
        else:
            marker = 'o'
        leg_label = ''
        for kplt in kplt_list:
            if j == Nrhoc-1:
                leg_label = 'K=%g, ' % Kentr_list[i] + labels[kplt]
            ax.scatter(prop[i, j, 5], prop[i, j, kplt], s=40, marker=marker,
                       color=color_list[i], label=leg_label)

# if kplt == 3:
#     ax.set_xscale('log')
#     ax.set_yscale('log')

ax.set_xlabel('rhoc')
ax.set_ylabel(labels[kplt])
ax.legend(loc=2, prop={'size': 25}, fancybox=True, framealpha=0.5)
pl.subplots_adjust(bottom=0.115, left=0.12, top=0.98, right=0.98)
fig.savefig(dir_main + labels[kplt] + '.png', dpi=300)

