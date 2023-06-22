import os
from funcs_tidal_eq import *
from para_tidal_eq import *

rhoname_initial = savedir + 'polytrope_profile_npoly%.5f.txt' % npoly
potname = savedir + 'potential.txt'
# assuming no stellar rotation and only considering quadrupolar tides

Niter = 0   # next iteration to be done (make sure potential%d and rho%d exist)
Nx, Ny, Nz, xarr, yarr, zarr = set_grid(Lmax, Nresz)

if Niter == 0:   # start from scratch
    # initial density profile
    rhoarr = map_LaneEmden(rhoname_initial, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhoc, Kentr)
    # rhoarr = map_LaneEmden_x0(rhoname_initial, Nx, Ny, Nz, xarr, yarr, zarr, npoly, x0surf, Kentr)
    # print('max(rho)=', np.amax(rhoarr))
else:
    rhoname = savedir + 'rho%d.txt' % Niter
    rhoarr = read_rho(rhoname, Nx, Ny, Nz)

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

# exit()

# update the density profile (including tidal potential)
# rhoarr = update_rho(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhoc, Kentr, Qbh, qstar)  # old coordinates
rhoarr = update_rho_sma(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, rhoc, Kentr, Qbh, qstar, sma)
# rhoarr = update_rho_x0(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, x0surf, Kentr, Qbh, qstar)  # bad!!
# print('max(rho)=', np.amax(rhoarr))
Niter += 1
write_rho(rhoarr, Niter, Nx, Ny, Nz)

qstar_old = qstar
qstar, xcom = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
frac_delta_q = (qstar - qstar_old)/qstar
print('qstar(%d)=%.5f, frac_delta_q=%.3e' % (Niter, qstar, frac_delta_q))

# exit()

# ----- check convergence based on total stellar mass
drifting = False
while abs(qstar - qstar_old)/qstar > rtol:
    os.system(srcdir + 'run %d %.2f %.2f %d' % (Niter, npoly, Lmax, Nresz))
    print('finished iteration %d' % Niter)
    Phiarr = read_Phi(potname, Nx, Ny, Nz)
    os.system('mv ' + potname + ' ' + savedir + 'potential%d.txt' % Niter)
    # rhoarr = update_rho(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhoc, Kentr, Qbh, qstar)
    rhoarr = update_rho_sma(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, rhoc, Kentr, Qbh, qstar, sma)
    # rhoarr = update_rho_x0(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, x0surf, Kentr, Qbh, qstar)  # very bad!
    Niter += 1
    write_rho(rhoarr, Niter, Nx, Ny, Nz)
    qstar_old, frac_delta_q_old = qstar, frac_delta_q
    qstar, xcom = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
    frac_delta_q = (qstar - qstar_old)/qstar
    print('qstar(%d)=%.5f, frac_delta_q=%.3e' % (Niter, qstar, frac_delta_q))
    # if Niter > 2 and abs(frac_delta_q) > abs(frac_delta_q_old):
    #     drifting = True
    #     break
# if drifting:
#     print('drifting has started, please use potential%d and rho%d as the final result' % (Niter-1, Niter-1))
# else:
#     print('solution has converged well because Mstar(%d)=Mstar(%d)!' % (Niter, Niter-1))
print('solution has converged well because Mstar(%d)=Mstar(%d)!' % (Niter, Niter-1))

# exit()

# calculate the position of the L1 point
# potential cloest to the x-axis (for iy = 0, iz = 0, but not exactly on x-axis)
Phixarr = [Phiarr[i, 0, 0] + Phitidal_sma(xarr[i], yarr[0], zarr[0], Qbh, qstar, sma) for i in range(Nx)]
xL1, xL2, PhiL1, PhiL2 = findL1L2(Phixarr, Nx, xarr)

# determine the position of stellar surface along x-axis
ixc = np.searchsorted(xarr, 0.)
isurf = ixc + 1
eps_small = 1e-10   # a small number equivalent to 0
while rhoarr[isurf, 0, 0] > eps_small:
    isurf += 1
xsurf = xarr[isurf-1]  # last point of non-zero density
print('****FINAL**** xL1=%.3f, PhiL1=%.3f, xsurf=%.3f' % (xL1, PhiL1, xsurf))

iL1 = np.searchsorted(xarr, xL1)
# note: xarr[iL1-1] < xL1 <= xarr[iL1]
if xsurf < xarr[iL1-1]:
    print('equilibrium star is detached!')
else:
    print('equilibrium star is filling up Roche Lobe!')
