import numpy as np
import os
from scipy.interpolate import interp1d
from math import pi
from multiprocessing import Process

# run a large number of equilibrium solutions for different Kentr and rhoc (fixing Qbh, sma, npoly)

# note: for Ncpu=12 on my desktop, it takes ~10 minutes for each Kentr

# ---- manually set
# Kentr_list = [0.1, 0.15, 0.2, 0.3]
# rhocmin, rhocmax = 5, 12

# ---- consider a particular entropy given by a low-mass main-sequence star
npoly = 1.5
ximax, phimax = 3.65375, 2.71409   # obtained from LaneEmden.py
Mstar = 0.5  # Msun
rhocmin, rhocmax = 2, 15
Rstar_fname = './data_figs/mass_radius_3Gyr.txt'
data = np.loadtxt(Rstar_fname, skiprows=1)
Marr, Rarr = data[:, 0], data[:, 1]
RM_interp = interp1d(Marr, Rarr)
Rstar = RM_interp(Mstar)   # Rsun
Kentr = (4*pi)**(1./npoly)/(npoly+1)*(Mstar/phimax)**(1-1./npoly)*(Rstar/ximax)**(3./npoly-1)
# print('Mstar, Rstar, Kentr = ', Mstar, Rstar, Kentr)
Kentr_list = [Kentr]

Ncpu = 12   # number of processors used (<=16 for best performance on desktop)
Nrhoc = 2*Ncpu
rhocarr = np.linspace(rhocmin, rhocmax, Nrhoc, endpoint=True)

# note: for large rhoc >~ 50, the star is small,
# so we reduce Lmax from the default value of 2.2 to 1.5


def run_tidal_eq(jlist, Kentr, s):
    # jlist is a chunk of range(Nrhoc), i for Kentr_list index, s is a random seed (not used)
    np.random.seed(s)
    for j in jlist:
        os.system('python ./tidal_equilibrium.py %.5f %.5f' % (Kentr, rhocarr[j]))


NK = len(Kentr_list)
jlist_chunks = np.array_split(range(Nrhoc), Ncpu)

if __name__ == '__main__':
    for i in range(NK):
        Kentr = Kentr_list[i]
        print('working on Kentr=', Kentr)
        Kentr_dir = './data_figs/Kentr%.3f' % Kentr
        if not os.path.exists(Kentr_dir):
            os.system('mkdir -p ' + Kentr_dir)
        for j in range(Nrhoc):
            rhoc = rhocarr[j]
            rhoc_dir = './data_figs/Kentr%.3f/rhoc%.3f' % (Kentr, rhoc)
            if os.path.exists(rhoc_dir):
                os.system('rm -rf ' + rhoc_dir)
            os.system('mkdir -p ' + rhoc_dir)
        procs = [Process(target=run_tidal_eq,
                         args=(jlist_chunks[n], Kentr, np.random.randint(10)))
                 for n in range(Ncpu)]
        for p in procs:
            p.start()
        for p in procs:
            p.join()
