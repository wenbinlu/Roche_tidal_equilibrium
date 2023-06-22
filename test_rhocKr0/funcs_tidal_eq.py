import numpy as np
from math import pi, floor, ceil, sqrt
from dir_info import *
from my_func import *


def set_grid(Lmax, Nresz):
    xarr = np.linspace(-Lmax, Lmax, 2*Nresz, endpoint=False)
    dx = xarr[1] - xarr[0]
    xarr += dx/2   # middle of the grid (not the same as used in Poisson3D.c)
    zarr = np.linspace(0, Lmax, Nresz, endpoint=False)
    dz = zarr[1] - zarr[0]
    zarr += dz/2   # same as used in Poisson3D.c
    return 2*Nresz, Nresz, Nresz, xarr, zarr, zarr


def map_LaneEmden(rhoname_initial, Nx, Ny, Nz, xarr, yarr, zarr, npoly, rhoc, Kentr):
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
    alpha_unit = sqrt((npoly+1)*Kentr*rhoc**(1./npoly - 1)/(4*pi))
    rmax = ximax*alpha_unit
    print('for initial profile, rmax=%.3f' % rmax)
    for i in range(Nx):
        x = xarr[i]
        for j in range(Ny):
            y = yarr[j]
            for k in range(Nz):
                z = zarr[k]
                r = sqrt(x*x + y*y + z*z)
                xi = r / alpha_unit
                # linear interpolation to obtain the density
                if xi >= ximax:
                    rhoarr_new[i, j, k] = 0.0
                    continue
                ixi = np.searchsorted(xiarr, xi)
                # note: xiarr[ixi-1] < xi < xiarr[ixi]
                slope = (thearr[ixi] - thearr[ixi-1])/(xiarr[ixi] - xiarr[ixi-1])
                the = thearr[ixi] + (xi - xiarr[ixi])*slope
                rhoarr_new[i, j, k] = the**npoly * rhoc
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


def Phitidal_sma(x, y, z, Qbh, qstar, sma):  # for a given SMA
    mu = Qbh / (Qbh + qstar)
    return -Qbh/sqrt((x-sma)**2 + y**2 + z**2) - (Qbh+qstar)/(2*sma**3) * ((x-mu*sma)**2 + y**2)


def update_rho_rhocKr0(Phiarr, Nx, Ny, Nz, xarr, yarr, zarr, xcom, npoly, rhoc, Kentr, x0surf, Qbh, qstar, atol):
    rhoarr_new = np.zeros((Nx, Ny, Nz), dtype=float)
    ixc = np.searchsorted(xarr, 0.)
    ix0 = np.searchsorted(xarr, x0surf)
    # note: xarr[ixc-1] < 0 <= xarr[ixc]
    xc, yc, zc = xarr[ixc], yarr[0], zarr[0]   # position very close to the center
    Phic = Phiarr[ixc, 0, 0]
    Phi0 = Phiarr[ix0, 0, 0]
    sma = findsma(Phic, Phi0, xc, yc, zc, npoly, rhoc, Kentr, x0surf, Qbh, qstar, atol)
    # set the surface of the star
    x0 = xarr[ix0]
    # rmax = x0   # maximum radial extension of the star
    Phitot0 = Phi0 + Phitidal_sma(x0, yc, zc, Qbh, qstar, sma)
    Phixarr = [Phiarr[i, 0, 0] + Phitidal_sma(xarr[i], yc, zc, Qbh, qstar, sma)
               for i in range(Nx)]
    xL1, PhiL1 = findL1(Phixarr, Nx, xarr)
    rmax = xL1 - xcom
    for k in range(Nz):
        z = zarr[k]
        if z >= rmax:
            break
        for j in range(Ny):
            y = yarr[j]
            if y >= rmax:
                break
            for i in range(Nx):
                x = xarr[i]
                if sqrt(x**2 + y**2 + z**2) > rmax:
                    continue
                Phitot = Phiarr[i, j, k] + Phitidal_sma(x, y, z, Qbh, qstar, sma)
                if Phitot > Phitot0:
                    continue
                rhoarr_new[i, j, k] = ((Phitot0 - Phitot)/(Kentr*(npoly+1)))**npoly
    return rhoarr_new


def stellar_mass(rhoarr, Nx, Ny, Nz, xarr, yarr, zarr):
    qstar, xqstar = 0., 0.
    dx = xarr[1] - xarr[0]
    dy = yarr[1] - yarr[0]
    dz = zarr[1] - zarr[0]
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                xqstar += xarr[i]*rhoarr[i, j, k] * dx*dy*dz
                qstar += rhoarr[i, j, k] * dx*dy*dz
    return qstar, xqstar/qstar


def findL1(Phixarr, Nx, xarr):
    ixc = np.searchsorted(xarr, 0.)
    # note: xarr[ixc-1] < 0 <= xarr[ixc]
    diffPhix = np.diff(Phixarr)
    iL1 = ixc+1
    while iL1 < Nx-2:
        if diffPhix[iL1] * diffPhix[iL1+1] < 0:
            break
        iL1 += 1
    if iL1 == Nx-2:
        print('no L1 point detected, please use a larger Lmax!')
        exit()
    iL1 += 1   # going back to the index of the xarr grid
    # print(xarr[iL1-1], xarr[iL1], xarr[iL1+1])
    # print(Phixarr[iL1-1], Phixarr[iL1], Phixarr[iL1+1])
    i = iL1     # now L1 point should be between iL1-1 and iL1+1
    Amatrix = np.array([[xarr[i-1]**2, xarr[i-1], 1.],
                        [xarr[i]**2, xarr[i], 1.],
                        [xarr[i+1] ** 2, xarr[i+1], 1.]])
    Amatrix_inv = np.linalg.inv(Amatrix)
    abc = np.dot(Amatrix_inv, [Phixarr[i-1], Phixarr[i], Phixarr[i+1]])
    xL1 = -abc[1]/(2*abc[0])
    PhiL1 = abc[2] - abc[1]**2/(4*abc[0])
    # print(xL1, PhiL1)
    return xL1, PhiL1


def funcsma(sma, lhs, xc, yc, zc, x0surf, Qbh, qstar):
    return Phitidal_sma(x0surf, yc, zc, Qbh, qstar, sma) - Phitidal_sma(xc, yc, zc, Qbh, qstar, sma) - lhs


def findsma(Phic, Phi0, xc, yc, zc, npoly, rhoc, Kentr, x0surf, Qbh, qstar, atol):
    lhs = rhoc**(1./npoly) * Kentr * (npoly + 1) + Phic - Phi0
    mu = Qbh/(Qbh + qstar)
    smaguess = abs(Qbh*x0surf**2*(1+0.5/mu)/lhs)**(1./3)   # based on quadrupolar tide only
    sma = bisec(funcsma, 0.5*smaguess, 1.5*smaguess, atol, lhs, xc, yc, zc, x0surf, Qbh, qstar)
    return sma


