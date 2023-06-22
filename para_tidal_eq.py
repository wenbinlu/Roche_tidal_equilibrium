# free parameters for "tidal_equilibrium.py"
rhoc = 14.  # need tests to see if it fills Roche Lobe
Kentr = 0.5    # fixed entropy

Qbh = 1e6
sma = 100.    # SMA
npoly = 1.5   # polytropic index

# accuracy of BC requires large Lmax, but grid resolution is better at small Lmax
Lmax = 2.2   # maximum box size (reasonable size: ~2.0)

Nresz = 100   # resolution per dimension (reasonable resolution: 100)
rtol = 1e-4   # fractional total-mass change (reasonable convergence: 1e-4)

# note1: frac mass difference is reduced by a factor of ~2 in each iteration
# note2: For some reason, once the solution has reasonably converged,
# it starts to drift away slowly. Currently, the convergence criterion detects
# the drifting behavior and stops the iteration once drifting starts
# Under this procedure, it typically takes <= 10 iterations to "converge"
# note3: the convergence is not very accurate (potentially because the finite box size/resolution)
