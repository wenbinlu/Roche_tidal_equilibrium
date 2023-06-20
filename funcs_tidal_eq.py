import numpy as np
from math import pi, floor, ceil, sqrt
from dir_info import *


def set_grid(Lmax, Nres):
    xarr = np.linspace(0, Lmax, Nres, endpoint=False)
    dx = xarr[1] - xarr[0]
    xarr += dx/2   # middle of the grid (as used in Poisson3D.c)
    return Nres, Nres, Nres, xarr, xarr, xarr


def map_LaneEmden(rhoname_initial, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhocbar):
    rhoarr_new = np.zeros((Nx, Ny, Nz), dtype=float)
    # map the LaneEmden solution into 3D grid and save the density profile
    with open(rhoname_initial, 'r') as f:
        Nxi = int(f.readline().strip('\n'))
        xiarr = np.zeros(Nxi, dtype=float)
        thearr = np.zeros(Nxi, dtype=float)
        for i in range(Nxi):
            row = f.readline().strip('\n').split(' ')
            xiarr[i] = float(row[0])
            thearr[i] = float(row[1])
    # note that xiarr is monotonically increasing and xiarr[0] = 0.0
    ximax = xiarr[-1]   # surface radius of the star
    rmax = ximax/sqrt(4*pi/rhocbar**(1./npoly-1))
    print('for initial profile, rmax=%.3f' % rmax)
    for i in range(Nx):
        x = xarr[i]
        for j in range(Ny):
            y = yarr[j]
            for k in range(Nz):
                z = zarr[k]
                r = sqrt(x*x + y*y + z*z)
                xi = r * sqrt(4*pi/rhocbar**(1./npoly-1))
                # linear interpolation to obtain the density
                if xi >= ximax:
                    rhoarr_new[i, j, k] = 0.0
                    continue
                ixi = np.searchsorted(xiarr, xi)
                # note: xiarr[ixi-1] < xi < xiarr[ixi]
                slope = (thearr[ixi] - thearr[ixi-1])/(xiarr[ixi] - xiarr[ixi-1])
                the = thearr[ixi] + (xi - xiarr[ixi])*slope
                rhoarr_new[i, j, k] = the**npoly * rhocbar
    return rhoarr_new


def write_rho(rhoarr_new, Niter, Nx, Ny, Nz):
    # write the density profile into a file to be used by C code
    with open(savedir + 'rho%d.txt' % Niter, 'w') as f:
        for i in range(Nx):
            for j in range(Ny):
                if (i+j) >= 1:
                    f.write('\n')
                for k in range(Nz):
                    if k >= 1:
                        f.write('\t')
                    f.write('%.8e' % rhoarr_new[i, j, k])


# read 3D density solution from a file
def read_rho(rhoname, Nx, Ny, Nz):
    rhoarr_new = np.zeros((Nx, Ny, Nz), dtype=float)
    with open(rhoname, 'r') as f:
        for i in range(Nx):
            for j in range(Ny):
                row = f.readline().strip('\n')
                if row == '':  # end of file
                    break
                dat = row.split('\t')
                rhoarr_new[i, j, :] = [float(rho) for rho in dat]
                # try:
                #     rhoarr_new[i, j, :] = [float(rho) for rho in dat]
                # except ValueError:
                #     print(len(dat), dat)
    return rhoarr_new


# read the potential solution from a file
def read_Phi(potname, Nx, Ny, Nz):
    Phiarr_new = np.zeros((Nx, Ny, Nz), dtype=float)
    with open(potname, 'r') as f:
        row = f.readline()
        while row != '':  # til the end
            if row[0] == 'j':
                j = int(row[4:].strip('\n').replace(' ', '')) - 1   # current y-index
                # print('j=%d' % j)
                k = 0   # counter for z-index
            if len(row) > 40:  # contains data
                dat = row.split('\t')
                # print(dat)
                # the last element of dat is '\n' (is removed)
                Phiarr_new[:, j, k] = [float(Phi) for Phi in dat[:-1]]
                k += 1
            row = f.readline()
    return Phiarr_new


def update_rho(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhocbar, lam2):
    rhoarr_new = np.zeros((Nx, Ny, Nz), dtype=float)
    Phic = Phiarr[0, 0, 0]
    for k in range(Nz):
        z = zarr[k]
        rhs_zaxis = rhocbar**(1./npoly) - (Phiarr[0, 0, k] - Phic + (1-lam2)*z**2)
        if rhs_zaxis <= 0:  # max surface extension on the z-axis
            break
        for j in range(Ny):
            y = yarr[j]
            rhs_yaxis = rhocbar**(1./npoly) - (Phiarr[0, j, k] - Phic
                                               + lam2 * y**2 + (1-lam2)*z**2)
            if rhs_yaxis <= 0:  # max surface extension on the y-axis (for given z)
                break
            for i in range(Nx):
                x = xarr[i]
                rhs = rhocbar**(1./npoly) - (Phiarr[i, j, k] - Phic
                                             -x**2 + lam2 * y**2 + (1-lam2) * z**2)
                if rhs <= 0:  # max extend of surface in the x-direction (for given z, y)
                    break
                rhonew = rhs**npoly
                if i >= 1 and rhonew > rhoarr_new[i-1, j, k]:
                    break  # do not allow increasing density along x-direction
                rhoarr_new[i, j, k] = rhonew
    # --- the following does not capture the cusp near L1 point
    # for i in range(Nx):
    #     x = xarr[i]
    #     for j in range(Ny):
    #         y = yarr[j]
    #         for k in range(Nz):
    #             z = zarr[k]
    #             rhs = rhocbar**(1./npoly) - (Phiarr[i, j, k] - Phic
    #                                          -x**2 + lam2 * y**2 + (1-lam2) * z**2)
    #             # rhs = rhocbar**(1./npoly) - (Phiarr[i, j, k] - Phic)
    #             if rhs <= 0:
    #                 break   # outside stellar surface
    #             rhonew = rhs**npoly
    #             # --- the following does not properly detect stellar surface
    #             # if rhonew < rhocbar:
    #             #     rhoarr_new[i, j, k] = rhonew
    #             # else:
    #             #     rhoarr_new[i, j, k] = 0.
    #             # ---- the following works well when star is far from RLO
    #             beyond_surface = False
    #             if i >= 1 and rhonew > rhoarr_new[i-1, j, k]:
    #                 beyond_surface = True
    #             if j >= 1 and rhonew > rhoarr_new[i, j-1, k]:
    #                 beyond_surface = True
    #             if k >= 1 and rhonew > rhoarr_new[i, j, k-1]:
    #                 beyond_surface = True
    #             # ---- only allow monotonically decreasing density
    #             if beyond_surface:
    #                 rhoarr_new[i, j, k] = 0.
    #             else:
    #                 rhoarr_new[i, j, k] = rhonew
    return rhoarr_new


def stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr):
    Mstar = 0.
    dx = xarr[1] - xarr[0]
    dy = yarr[1] - yarr[0]
    dz = zarr[1] - zarr[0]
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                Mstar += rhoarr[i, j, k] * dx*dy*dz
    return Mstar