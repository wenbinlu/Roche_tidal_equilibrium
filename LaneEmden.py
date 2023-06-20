import numpy as np
from dir_info import *

savename = 'polytrope_profile'

npoly = 1.5
res = 1e-5   # resolution in dxi, method is only first-order accurate

# boundary condition
xi = 0.0 + 1e-15
the = 1.0
phi = 0.0

# initialization
xiarr = [xi]
thearr = [the]
phiarr = [phi]

# resolution (using a first-order Euler scheme)
dxi = res
ximax, phimax = 0, 0.   # only initialization

while the > 0.0:
    dthedxi = -phi/xi**2
    dphidxi = xi**2 * the**npoly
    if the + dthedxi * dxi < 0:  # take a partial step
        ximax = xi - the/dthedxi
        dxi = ximax - xi
        phimax = phi + dphidxi * dxi
        break
    the += dthedxi * dxi
    phi += dphidxi * dxi
    xi += dxi

    xiarr += [xi]
    thearr += [the]
    phiarr += [phi]
    

print('ximax=%.5f, phimax=%.5f' % (ximax, phimax))

# save a low-resolution xiarr, thearr, phiarr profiles in a data file
Nxi = 1000
savexiarr = np.linspace(0, xiarr[-1], Nxi, endpoint=False)
savethearr = np.zeros(Nxi, dtype=float)
savephiarr = np.zeros(Nxi, dtype=float)
j = 0
for i in range(len(xiarr)):
    xi = xiarr[i]
    if xi >= savexiarr[j]:
        the_slope = (thearr[i+1] - thearr[i]) / (xiarr[i+1] - xi)
        phi_slope = (phiarr[i+1] - phiarr[i]) / (xiarr[i+1] - xi)
        the = thearr[i] + the_slope * (savexiarr[j] - xi)
        phi = phiarr[i] + phi_slope * (savexiarr[j] - xi)
        savethearr[j] = the
        savephiarr[j] = phi
        j += 1
        if j == Nxi:
            break
with open(savedir + savename + '.txt', 'w') as f:
    f.write('%08d' % len(savexiarr))
    for i in range(Nxi):
        f.write('\n%.10e %.10e %.10e' % (savexiarr[i], savethearr[i], savephiarr[i]))
