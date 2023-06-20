# free parameters for "tidal_equilibrium.py"
rhocbar = 9.9  # need tests to see if it fills Roche Lobe
lam2 = 0.3   # must be in the range (0, 1)

npoly = 1.5   # polytropic index
# accuracy of BC requires large Lmax, but grid resolution is better at small Lmax
Lmax = 2.2   # maximum box size (reasonable size: ~2.0)

Nres = 100   # resolution per dimension (reasonable resolution: 100)
tol = 1e-4   # fractional total-mass change (reasonable convergence: 1e-4)
# note: frac mass difference is reduced by a factor of ~2 in each iteration
# it typically takes <= 10 iterations to coverge
