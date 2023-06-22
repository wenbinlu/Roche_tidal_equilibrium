# free parameters for "tidal_equilibrium.py"
rhoc = 10.0  # need tests to see if it fills Roche Lobe
x0surf = 0.8  # fix one surface position along x axis (make sure |x0| < xL1)
Kentr = 0.5    # fixed entropy

Qbh = 1.0e3
# sma = 90.0    # SMA
npoly = 1.5   # polytropic index

# accuracy of BC requires large Lmax, but grid resolution is better at small Lmax
Lmax = 3.0   # maximum box size (reasonable size: ~2.5)

Nresz = 100   # resolution per dimension (reasonable resolution: 100)
tol = 1e-5   # fractional total-mass change (reasonable convergence: 1e-4)
atol = 1e-8   # absolute tolerance for bisection method (used in findsma)
# note: frac mass difference is reduced by a factor of ~2 in each iteration
# it typically takes <= 10 iterations to converge
