import matplotlib.pyplot as pl
from my_func import pltimg, round_sig
from dir_info import *
from funcs_tidal_eq import *
from para_tidal_eq import *
import sys

plt_case = int(sys.argv[1])    # 0 for pot and 1 for rho
Niter = int(sys.argv[2])    # [0 to max_Niter] pick the density profile to plot

# for the potential plot, it is best to pick Niter = max_Niter (fyr self-consistency)

# ----- plot a slice at a given z
z_plt = 0.   # must be less than Lmax
Ncont = 15   # number of contours for potential plot only

# to get rid of zeros
eps_small = 1e-10

potname = savedir + 'potential.txt'
rhoname = savedir + 'rho%d.txt' % Niter

if plt_case == 0:
    savename = 'fig_pot'
else:
    savename = 'fig_rho%d' % Niter

rhoCB_levels = [-4, -3, -2, -1, 0, 0.5, 1]
while rhoCB_levels[-1] > np.log10(rhocbar):
    del rhoCB_levels[-1]
rhoCB_ticklabels = [('%g' % num) for num in rhoCB_levels]

Nx, Ny, Nz, xarr, yarr, zarr = set_grid(Lmax, Nres)

# read the converged dimensionless potential \bar{Phi}
Phiarr = read_Phi(potname, Nx, Ny, Nz)
rhoarr = read_rho(rhoname, Nx, Ny, Nz)

# get rid of all the zeros in density profile
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            rhoarr[i, j, k] = max(eps_small, rhoarr[i, j, k])

k_plt = min(Nz-1, max(0, np.searchsorted(zarr, z_plt)-1))
# such that zarr[k_plt-1] < z_plt < zarr[k_plt]
# ---- choose what to plot
if plt_case == 0:
    # add the tidal potential
    for i in range(Nx):
        x = xarr[i]
        for j in range(Ny):
            y = yarr[j]
            for k in range(Nz):
                z = zarr[k]
                Phitidal = -x**2 + lam2 * y**2 + (1-lam2) * z**2
                Phiarr[i, j, k] += Phitidal
    pltarr = Phiarr[:, :, k_plt]
    pltlabel = r'$\Phi(z=%.1f)$' % zarr[k_plt]
    if k_plt == 0:  # in the x-y plane
        # find the L1 point
        j, k = 0, 0
        diffPhi = np.diff(Phiarr[:, j, k])
        iL1 = 0
        while iL1 < Nx-1:
            if diffPhi[iL1] * diffPhi[iL1+1] < 0:
                break
            iL1 += 1
        xL1 = xarr[iL1]
        # find the stellar surface along x-axis
        isurf = 0
        while rhoarr[isurf, 0, 0] > eps_small:
            isurf += 1
        xsurf = xarr[isurf-1]
else:
    pltarr = np.log10(rhoarr[:, :, k_plt])
    pltlabel = r'$\log\rho$'

xlabel = r'$x$'
ylabel = r'$y$'

fig = pl.figure(figsize=(13, 10))
ax = fig.add_axes([0.12, 0.10, 0.87, 0.87])

min_val, max_val = np.amin(pltarr), np.amax(pltarr)
print('min, max: ', min_val, max_val)
step = round_sig((max_val-min_val)/Ncont, sig=1)  # round to sig figure = 1
res = step
potCB_levels = np.arange(ceil(min_val/res)*res,
                      floor(max_val/res)*res + step, step)
potCB_ticklabels = [('%.1f' % num) for num in potCB_levels]

if plt_case == 0:
    pltimg(ax, xarr, yarr, pltarr, xlabel, ylabel, pltlabel, 'bwr',
           potCB_levels, potCB_ticklabels, flag_contour=True)
    if k_plt == 0:  # show the L1 point
        ax.plot(xL1, yarr[0], 'ko', ms=5, alpha=1, fillstyle='none')  # empty
        ax.plot(xsurf, yarr[0], 'ko', ms=5, alpha=1)  # filled symbol
    # overplot the density contours
    X, Y = np.meshgrid(xarr, yarr)
    logrhoarr = np.log10(rhoarr[:, :, k_plt])
    CS = ax.contour(X, Y, logrhoarr.transpose(),
                    rhoCB_levels, linestyles='solid',
                    colors='orange', linewidths=2, alpha=1)
    fmt = {}
    for l, s in zip(CS.levels, rhoCB_ticklabels):
        fmt[l] = s
    pl.clabel(CS, CS.levels, inline=True, fmt=fmt,
              fontsize=30, colors=None)
else:
    pltimg(ax, xarr, yarr, pltarr, xlabel, ylabel, pltlabel, 'bwr',
           rhoCB_levels, rhoCB_ticklabels, flag_contour=True)

# pl.subplots_adjust(bottom=0.13, left=0.12, top=0.98, right=0.98)
pl.savefig(savedir + savename + '.png')
