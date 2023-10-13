import numpy as np
import os
from scipy.interpolate import interp1d
from math import pi
from multiprocessing import Process
from para_tidal_eq import *
from dir_info import *

# run a large number of equilibrium solutions for different Kentr and rhoc
# (fixing Qbh, sma, npoly, as specified by para_tidal_eq.py)

# note: for Ncpu=12 on my desktop, it takes ~10 minutes for each Kentr

rhocmin, rhocmax = 1.0, 5.0

# ---- manually set a list of entropies
# Kentr_list = [0.1, 0.15, 0.2, 0.3]

# ---- consider a particular entropy given by a star of given mass and radius
Mstar = 2.0  # Msun
Rstar = 2.2

# ---- use main-sequence star mass-radius relation
# Rstar_fname = savedir + 'mass_radius_3Gyr.txt'
# data = np.loadtxt(Rstar_fname, skiprows=1)
# Marr, Rarr = data[:, 0], data[:, 1]
# RM_interp = interp1d(Marr, Rarr)
# Rstar = RM_interp(Mstar)   # Rsun

LaneEmden_fname = 'polytrope_profile_npoly%.5f' % npoly + '.txt'
if not os.path.exists(savedir + LaneEmden_fname):
    # need to run LaneEmden.py to create the polytrope_profile
    os.system('python ' + pydir + 'LaneEmden.py')
data = np.loadtxt(savedir + LaneEmden_fname, skiprows=1)
ximax, phimax = data[-1, 0], data[-1, 2]

Kentr = (4*pi)**(1./npoly)/(npoly+1)*(Mstar/phimax)**(1-1./npoly)*(Rstar/ximax)**(3./npoly-1)
# print('Mstar, Rstar, Kentr = ', Mstar, Rstar, Kentr)
# exit()

Kentr_list = [Kentr]
Ncpu = 12   # number of processors used (<=16 for best performance on desktop)
Nrhoc = 2*Ncpu   # each processor calculates two cases
rhocarr = np.linspace(rhocmin, rhocmax, Nrhoc, endpoint=True)

# note: for very large rhoc, the star is small,
# so we need to reduce Lmax from the default value of 2.0 to <~1


def run_tidal_eq(jlist, Kentr, s):
    # jlist is a chunk of range(Nrhoc), i for Kentr_list index, s is a random seed (not used)
    np.random.seed(s)
    for j in jlist:
        os.system('python ' + dir_main + 'tidal_equilibrium.py %.5f %.5f' % (Kentr, rhocarr[j]))


NK = len(Kentr_list)
jlist_chunks = np.array_split(range(Nrhoc), Ncpu)

# creat all the relevant directories
for i in range(NK):
    Kentr = Kentr_list[i]
    Kentr_dir = savedir + 'Kentr%.3f' % Kentr
    if not os.path.exists(Kentr_dir):
        os.system('mkdir -p ' + Kentr_dir)
    for j in range(Nrhoc):
        rhoc = rhocarr[j]
        rhoc_dir = savedir + 'Kentr%.3f/rhoc%.3f' % (Kentr, rhoc)
        if os.path.exists(rhoc_dir):
            os.system('rm -rf ' + rhoc_dir)
        os.system('mkdir -p ' + rhoc_dir)


if __name__ == '__main__':
    for i in range(NK):
        Kentr = Kentr_list[i]
        print('working on Kentr=', Kentr)
        procs = [Process(target=run_tidal_eq,
                         args=(jlist_chunks[n], Kentr, np.random.randint(10)))
                 for n in range(Ncpu)]
        for p in procs:
            p.start()
        for p in procs:
            p.join()
