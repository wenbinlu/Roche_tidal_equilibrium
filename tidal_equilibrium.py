import os
from funcs_tidal_eq import *
from para_tidal_eq import *

rhoname_initial = savedir + 'polytrope_profile.txt'
potname = savedir + 'potential.txt'
# assuming no stellar rotation and only considering quadrupolar tides

Niter = 0   # number of iterations so far
Nx, Ny, Nz, xarr, yarr, zarr = set_grid(Lmax, Nres)

# initial dimensionless density \bar{rho}
rhoarr = map_LaneEmden(rhoname_initial, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhocbar)

Mstar = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
print('Mstar(%d)=%.5f' % (Niter, Mstar))
write_rho(rhoarr, Niter, Nx, Ny, Nz)

# compile the C code
compile_command = 'gcc -o run main.c boundary_functions.c parameters_and_boundary_types.c ../lib/libutil.a'
os.chdir(srcdir)
os.system(compile_command)

# run the C code
os.system(srcdir + 'run %d %.2f %.2f %d' % (Niter, npoly, Lmax, Nres))
print('finished iteration %d' % Niter)
# read dimensionless potential \bar{Phi}
Phiarr = read_Phi(potname, Nx, Ny, Nz)

# update the density profile (including tidal potential)
rhoarr = update_rho(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhocbar, lam2)
Niter += 1
write_rho(rhoarr, Niter, Nx, Ny, Nz)

Mstar_old = Mstar
Mstar = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
print('Mstar(%d)=%.5f, frac_Delta_M=%.3e' % (Niter, Mstar, (Mstar - Mstar_old)/Mstar))

# ----- check convergence based on total stellar mass
while abs(Mstar - Mstar_old)/Mstar > tol:
    os.system(srcdir + 'run %d %.2f %.2f %d' % (Niter, npoly, Lmax, Nres))
    print('finished iteration %d' % Niter)
    Phiarr = read_Phi(potname, Nx, Ny, Nz)
    rhoarr = update_rho(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhocbar, lam2)
    Niter += 1
    Mstar_old = Mstar
    Mstar = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
    print('Mstar(%d)=%.5f, frac_Delta_M=%.3e' % (Niter, Mstar, (Mstar - Mstar_old)/Mstar))
    write_rho(rhoarr, Niter, Nx, Ny, Nz)

print('solution has converged because Mstar(%d)=Mstar(%d)!' % (Niter, Niter-1))

# calculate the position of the L1 point
Phitot = np.zeros_like(xarr)   # only along the x-axis (j=k=0, not exactly y=z=0)
for i in range(Nx):
    x = xarr[i]
    Phitidal = -x**2 + lam2 * yarr[0]**2 + (1-lam2) * zarr[0]**2
    Phitot[i] = Phiarr[i, 0, 0] + Phitidal
diffPhi = np.diff(Phitot)   # gradient along the x-axis
iL1 = 0
while iL1 <= Nx-2:
    if diffPhi[iL1] * diffPhi[iL1+1] < 0:
        break
    iL1 += 1
# we find that diffPhi[iL1] < 0 and diff[iL1+1] > 0,
# where diffPhi[iL1] = Phitot[iL1+1] - Phitot[iL1]
# thus, the point closest to the potential minimum is x = xarr[iL1]
xL1 = xarr[iL1]

# determine the position of stellar surface along x-axis
isurf = 0
eps_small = 1e-10   # a small number equivalent to 0
while rhoarr[isurf, 0, 0] > eps_small:
    isurf += 1
xsurf = xarr[isurf-1]  # last point of non-zero density

if xsurf < xL1:
    print('equilibrium star is detached from Roche Lobe!')
else:
    print('equilibrium star is filling up Roche Lobe!')
