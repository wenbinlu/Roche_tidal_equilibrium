import os
from funcs_tidal_eq import *
from para_tidal_eq import *

rhoname_initial = savedir + 'polytrope_profile_npoly%.5f.txt' % npoly
potname = savedir + '../../data_figs/potential.txt'
# assuming no stellar rotation and only considering quadrupolar tides

Niter = 0   # number of iterations so far
Nx, Ny, Nz, xarr, yarr, zarr = set_grid(Lmax, Nresz)

# initial density profile
rhoarr = map_LaneEmden(rhoname_initial, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhoc, Kentr)
# print('max(rho)=', np.amax(rhoarr))

qstar, xcom = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
print('qstar(%d)=%.5f' % (Niter, qstar))
write_rho(rhoarr, Niter, Nx, Ny, Nz)

# compile the C code
compile_command = 'gcc -o run main.c boundary_functions.c parameters_and_boundary_types.c ../lib/libutil.a'
os.chdir(srcdir)
os.system(compile_command)
# print(compile_command)

# run the C code
run_command = srcdir + 'run %d %.2f %.2f %d' % (Niter, npoly, Lmax, Nresz)
os.system(run_command)
# print(run_command)
print('finished iteration %d' % Niter)
# read dimensionless potential \bar{Phi}
Phiarr = read_Phi(potname, Nx, Ny, Nz)
os.system('mv ' + potname + ' ' + savedir + 'potential%d.txt' % Niter)

# update the density profile (including tidal potential)
rhoarr = update_rho_rhocKr0(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, rhoc, Kentr, x0surf, Qbh, qstar, atol)
# print('max(rho)=', np.amax(rhoarr))
Niter += 1
write_rho(rhoarr, Niter, Nx, Ny, Nz)

qstar_old = qstar
qstar, xcom = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
print('Mstar(%d)=%.5f, frac_Delta_M=%.3e' % (Niter, qstar, (qstar - qstar_old)/qstar))


# ----- check convergence based on total stellar mass
while abs(qstar - qstar_old)/qstar > tol:
    os.system(srcdir + 'run %d %.2f %.2f %d' % (Niter, npoly, Lmax, Nresz))
    print('finished iteration %d' % Niter)
    Phiarr = read_Phi(potname, Nx, Ny, Nz)
    os.system('mv ' + potname + ' ' + savedir + 'potential%d.txt' % Niter)
    rhoarr = update_rho_rhocKr0(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, rhoc, Kentr, x0surf, Qbh, qstar, atol)
    Niter += 1
    qstar_old = qstar
    qstar, xcom = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
    print('Mstar(%d)=%.5f, frac_Delta_M=%.3e' % (Niter, qstar, (qstar - qstar_old)/qstar))
    write_rho(rhoarr, Niter, Nx, Ny, Nz)

print('solution has converged because Mstar(%d)=Mstar(%d)!' % (Niter, Niter-1))


# calculate the position of the L1 point
ixc = np.searchsorted(xarr, 0.)
# note: xarr[ixc-1] < 0 <= xarr[ixc]
xc, yc, zc = xarr[ixc], yarr[0], zarr[0]   # position very close to the center
mu = Qbh/(Qbh + qstar)
Phic = Phiarr[ixc, 0, 0]
ix0 = np.searchsorted(xarr, x0surf)
Phi0 = Phiarr[ix0, 0, 0]
sma = findsma(Phic, Phi0, xc, yc, zc, npoly, rhoc, Kentr, x0surf, Qbh, qstar, atol)
# potential clost to the x-axis (for iy = 0, iz = 0, but not exactly on x-axis)
Phixarr = [Phiarr[i, 0, 0] + Phitidal_sma(xarr[i], yc, zc, Qbh, qstar, sma) for i in range(Nx)]
xL1, PhiL1 = findL1(Phixarr, Nx, xarr)

# determine the position of stellar surface along x-axis
isurf = ixc + 1
eps_small = 1e-10   # a small number equivalent to 0
while rhoarr[isurf, 0, 0] > eps_small:
    isurf += 1
xsurf = xarr[isurf-1]  # last point of non-zero density
print('****FINAL**** xL1=%.3f, PhiL1=%.3f, xsurf=%.3f' % (xL1, PhiL1, xsurf))

iL1 = np.searchsorted(xarr, xL1)
# note: xarr[iL1-1] < xL1 <= xarr[iL1]
if xsurf < xarr[iL1-1]:
    print('equilibrium star is detached from Roche Lobe!')
else:
    print('equilibrium star is filling up Roche Lobe!')
