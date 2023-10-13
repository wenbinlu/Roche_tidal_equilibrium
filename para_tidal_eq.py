# rhoc and Kentr are specified by input arguments for tidal_equilibrium.py

# the following are fixed
Qbh = 1.4     # BH mass in Msun
sma = 6.0    # SMA in Rsun
npoly = 3.0   # polytropic index

# Qbh = 1.0e6     # BH mass in Msun
# sma = 300    # SMA in Rsun
# npoly = 1.5   # polytropic index

# --- accuracy of BC requires large Lmax, but grid resolution is better at small Lmax
Lmax = 4.0   # maximum box size for Poisson solver (reasonable size: ~2.0)

Nresz = 100   # resolution per dimension (reasonable resolution: 100 to 150)

# related to iterations
OnlySaveLast = True   # if true, only save the converged potential%d and rho%d files for each run
Niter_max = 15   # maximum number of iterations
rtol = 1e-4   # fractional total-mass change (reasonable convergence: 1e-4)

# note: fractional mass difference between adjacent iterations
# is reduced by a factor of ~2 in each iteration
