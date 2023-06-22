import numpy as np
from math import sqrt, sin, cos, asin, pi, acos, log10, floor, ceil, exp, log
import pylab as pl


def pltimg(ax, xarr, yarr, zarr, xlabel, ylabel, zlabel, cmap,
           CB_levels, CB_ticklabels, flag_contour):
    min_val, max_val = np.amin(zarr), np.amax(zarr)
    ax.set_xlabel(xlabel, labelpad=-2)
    ax.set_ylabel(ylabel)
    im = ax.imshow(zarr.transpose(),
                   interpolation='bicubic', origin='lower',
                   cmap=cmap, aspect='equal', alpha=0.7,
                   extent=(min(xarr), max(xarr),
                           min(yarr), max(yarr)))
    im.set_clim(vmin=min_val, vmax=max_val)

    CB = pl.colorbar(im, ax=ax, ticks=CB_levels, orientation='horizontal')
    CB.ax.set_xticklabels(CB_ticklabels)
    CB.ax.set_xlabel(zlabel, labelpad=3)
    # CB.ax.set_yticklabels(CB_ticklabels)
    # CB.ax.set_ylabel(zlabel, labelpad=3)
    CB.ax.minorticks_off()

    if flag_contour:
        X, Y = np.meshgrid(xarr, yarr)
        CS = ax.contour(X, Y, zarr.transpose(),
                        CB_levels, linestyles='solid',
                        colors='k', linewidths=2, alpha=0.5)
        fmt = {}
        for l, s in zip(CS.levels, CB_ticklabels):
            fmt[l] = s
        pl.clabel(CS, CS.levels, inline=True, fmt=fmt,
                  fontsize=30, colors=None)


def round_sig(x, sig=2):
    # round a number to a given significant digit
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


def bisec(y, xleft, xright, atol, *args):
    # use bisection method to find xleft<x<xright such that y(x) = 0
    # function y(x) must be monotonic (increasing or decreasing)
    yleft = y(xleft, *args)
    yright = y(xright, *args)
    # print(xleft, xright, *args)
    # return 0

    if yleft * yright > 0:
        print('bisection fails for parameters', xleft, xright, atol, *args)
        return 0
    while abs(xright-xleft) > atol:
        xmid = 0.5*(xright+xleft)
        ymid = y(xmid, *args)
        if ymid*yleft > 0:
            xleft = xmid
            yleft = ymid
        else:
            xright = xmid
    return 0.5*(xright+xleft)

