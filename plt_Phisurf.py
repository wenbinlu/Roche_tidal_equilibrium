import numpy as np
import pylab as pl
import os
from scipy.interpolate import interp1d
from math import pi, log
from glob import glob

Qbh, sma = 1e6, 300   # fixed parameters
# Qbh, sma = 1.0, 1.5

# these need to be the same as "run_Kentr_rhoc.py"
# Kentr_list = [0.2, 0.3, 0.4, 0.5]
# rhocmin, rhocmax = 10, 50
# Nrhoc = 24

qstar_target = 0.5

# Kentr_list = [0.1, 0.15, 0.2, 0.3]
Kentr_list = [0.157]
Nrhoc = 24

# ---- ximax and phimax only used to estimate Rstar for spherical star
npoly, ximax, phimax = 1.5, 3.65375, 2.71409   # obtained from LaneEmden.py

dir_main = '/Users/wenbinlu/PycharmProjects/Roche_tidal_equilibrium/data_figs/sma%.1f/' % sma
color_list = ['tomato', 'royalblue', 'olive', 'orange']
labels = ['qstar', 'xL1', 'PhitotL1+1.5Q/a', 'xsurf', 'Phitotsurf+1.5Q/a', 'rhoc']
# labels = ['qstar', 'xL1', 'PhitotL1', 'xsurf', 'Phitotsurf', 'rhoc']
kplt_list = [0]   # which property to be plotted [0-5]

NK = len(Kentr_list)
status = np.zeros((NK, Nrhoc), dtype=bool)
prop = np.zeros((NK, Nrhoc, 6), dtype=float)   # five properties of the equilibrium solution
maskrhoc = np.zeros((NK, Nrhoc), dtype=bool)   # mask unused simulations
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
            if 'converged' not in f.read():  # this simulation isn't useful
                maskrhoc[i, j] = True
                continue
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
    Kentr = Kentr_list[i]
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
        if maskrhoc[i, j]:
            continue  # skip the masked rhoc
        rhoc = prop[i, j, 5]
        qstar = prop[i, j, 0]
        Rstar = ((Kentr*(npoly+1))**npoly/(4*pi) *
                 (qstar/phimax)**(1.-npoly))**(1/(3.-npoly)) * ximax
        RLeff = 0.49*sma/(0.6+(Qbh/qstar)**(2./3)*log(1 + (Qbh/qstar)**(-1./3)))
        leg_label = ''
        # if j == Nrhoc-1:
        #     leg_label = r'$R_*/R_{\rm L}$'
        ax.scatter(rhoc, Rstar/RLeff, s=40, marker='*', color=color_list[i], label=leg_label)
        # ax.scatter(qstar, Rstar/RLeff, s=40, marker='*', color='k', label=leg_label)
        if status[i, j]:
            marker = 'x'
        else:
            marker = 'o'
        leg_label = ''
        for kplt in kplt_list:
            if j == 0:
                leg_label = labels[kplt] + ' (K=%g)' % Kentr_list[i]
            ax.scatter(rhoc, prop[i, j, kplt], s=40, marker=marker,
                       color=color_list[i], label=leg_label)

# if kplt == 3:
#     ax.set_xscale('log')
#     ax.set_yscale('log')
ax.grid()
ax.set_xlabel('rhoc')
# ax.set_ylabel(labels[kplt_list[-1]])
ax.legend(loc=2, prop={'size': 25}, fancybox=True, framealpha=0.5)
pl.subplots_adjust(bottom=0.115, left=0.12, top=0.98, right=0.98)
fig.savefig(dir_main + labels[kplt_list[-1]] + '.png', dpi=300)

