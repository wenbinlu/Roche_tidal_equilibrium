# the following are fixed
Qbh = 1e6     # BH mass in Msun
sma = 70.    # SMA in Rsun
npoly = 1.5   # polytropic index
# --- accuracy of BC requires large Lmax, but grid resolution is better at small Lmax
Lmax = 1.5   # maximum box size (reasonable size: ~2.0)

Nresz = 100   # resolution per dimension (reasonable resolution: 100)
rtol = 1e-5   # fractional total-mass change (reasonable convergence: 1e-4)

# note: frac mass difference is reduced by a factor of ~2 in each iteration
