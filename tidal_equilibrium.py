import os
import sys
from funcs_tidal_eq import *
from para_tidal_eq import *
from dir_info import *

# example run
# python tidal_equilibrium.py 0.578 2.043

# input parameters
Kentr = float(sys.argv[1])    #  # entropy constant [in units such that G=Msun=Rsun=1] 
rhoc = float(sys.argv[2])     # peak density [Msun/Rsun^3]  

Niter = 0   # start from this iteration number (if >= 1, make sure potential%d and rho%d exist)

savedir_Krhoc = savedir + 'Kentr%.3f/rhoc%.3f/' % (Kentr, rhoc)   # where data and prints are saved
if not os.path.exists(savedir_Krhoc):
    os.system('mkdir -p ' + savedir_Krhoc)   # '-p' allows a multi-layer directory to be created
log_file = open(savedir_Krhoc + "output.txt", "w")
sys.stdout = log_file
potname_C_output = savedir_Krhoc + 'potential.txt'   # keep this the same for all iterations

rhoname = savedir_Krhoc + 'rho%d.txt' % Niter
Nx, Ny, Nz, xarr, yarr, zarr = set_grid(Lmax, Nresz)

if Niter == 0:   # start from scratch
    # initial density profile
    rhoname_initial = dir_main + 'data_figs/polytrope_profile_npoly%.5f.txt' % npoly
    rhoarr = map_LaneEmden(rhoname_initial, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhoc, Kentr)
    # print('max(rho)=', np.amax(rhoarr))
    write_rho(rhoname, rhoarr, Nx, Ny, Nz)
else:
    rhoarr = read_rho(rhoname, Nx, Ny, Nz)

qstar, xcom = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
print('qstar(%d)=%.5f' % (Niter, qstar))

# compile the C code
compile_command = 'gcc -o run main.c boundary_functions.c parameters_and_boundary_types.c ../lib/libutil.a'
os.chdir(srcdir)
os.system(compile_command)

# run the C code
run_command = srcdir + 'run %d %.3f %.3f %d %.3f %.3f' % (Niter, npoly, Lmax, Nresz, rhoc, Kentr)
os.system(run_command)
# print(run_command)
print('finished iteration %d' % Niter)
# read dimensionless potential Phi
Phiarr = read_Phi(potname_C_output, Nx, Ny, Nz)
potname = savedir_Krhoc + 'potential%d.txt' % Niter
os.system('mv ' + potname_C_output + ' ' + potname)

# update the density profile (including tidal potential)
rhoarr = update_rho_sma(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, rhoc, Kentr, Qbh, qstar, sma)
# print('max(rho)=', np.amax(rhoarr))

qstar_old = qstar
qstar, xcom = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
frac_delta_q = (qstar - qstar_old)/qstar
print('qstar(%d)=%.5f, frac_delta_q=%.3e' % (Niter, qstar, frac_delta_q))

# ----- check convergence based on total stellar mass
while abs(qstar - qstar_old)/qstar > rtol:
    Niter += 1
    if Niter > Niter_max:
        break
    rhoname_old = rhoname
    if OnlySaveLast:  # remove rhoname_old file
        if os.path.exists(rhoname_old):
            os.system('rm ' + rhoname_old)
    rhoname = savedir_Krhoc + 'rho%d.txt' % Niter
    write_rho(rhoname, rhoarr, Nx, Ny, Nz)
    run_command = srcdir + 'run %d %.3f %.3f %d %.3f %.3f' % (Niter, npoly, Lmax, Nresz, rhoc, Kentr)
    os.system(run_command)
    print('finished iteration %d' % Niter)
    Phiarr = read_Phi(potname_C_output, Nx, Ny, Nz)  # always read the freshly generated C_output
    potname_old = potname
    if OnlySaveLast:  # remove files generated in intermediate interation steps
        if os.path.exists(potname_old):
            os.system('rm ' + potname_old)
    potname = savedir_Krhoc + 'potential%d.txt' % Niter
    os.system('mv ' + potname_C_output + ' ' + potname)
    rhoarr = update_rho_sma(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, rhoc, Kentr, Qbh, qstar, sma)
    # rhoarr = update_rho(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhoc, Kentr, Qbh, qstar)  # old units
    # rhoarr = update_rho_x0(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, x0surf, Kentr, Qbh, qstar)  # bad!
    qstar_old, frac_delta_q_old = qstar, frac_delta_q
    qstar, xcom = stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr)
    frac_delta_q = (qstar - qstar_old)/qstar
    print('qstar(%d)=%.5f, frac_delta_q=%.3e' % (Niter, qstar, frac_delta_q))

if Niter > Niter_max:
    print('maximum number of iterations (%d) reached, and the solution does not converge!' % Niter_max)
else:
    print('solution has converged because Mstar(%d)=Mstar(%d)!' % (Niter+1, Niter))

# exit()

# calculate the position of the L1 point
# potential closest to the x-axis (for iy = 0, iz = 0, but not exactly on x-axis)
Phixarr = [Phiarr[i, 0, 0] + Phitidal_sma(xarr[i], yarr[0], zarr[0], Qbh, qstar, sma) for i in range(Nx)]
xL1, xL2, PhiL1, PhiL2 = findL1L2(Phixarr, Nx, xarr)

# determine the position of stellar surface along x-axis
ixc = np.searchsorted(xarr, 0.)
isurf = ixc + 1
eps_small = 1e-10   # a small number equivalent to 0
while rhoarr[isurf, 0, 0] > eps_small:
    isurf += 1
xsurf = xarr[isurf-1]  # last point of non-zero density
Phitotsurf = Phixarr[isurf-1]
print('equilibrium result: qstar=%.5f, xL1=%.5f, PhitotL1+1.5Q/a=%.8f, xsurf=%.5f, Phitotsurf+1.5Q/a=%.8f, rhoc=%.5f'
      % (qstar, xL1, PhiL1, xsurf, Phitotsurf, rhoc))

iL1 = np.searchsorted(xarr, xL1)
# note: xarr[iL1-1] < xL1 <= xarr[iL1]
if xsurf < xarr[iL1-1]:
    print('FINAL: detached')
else:
    print('FINAL: overflowing')

# remove grid file (not used)
gridfile = savedir_Krhoc + 'GridX.txt'
if os.path.exists(gridfile):
    os.system('rm ' + gridfile)

log_file.close()
