# Roche_tidal_equilibrium
The goal is to calculate the density profile of a polytropic star in equilibrium under Newtonian tidal potential of a point-mass companion. The current implementation assumes a star in synchronous rotation and calculations are done in the corotating frame with the coordinate origin at the center of the star.

There are two major pieces of codes: a C-based Poisson solver (poisson3D) from [Derek Harrison's Github repository](https://github.com/derekharrison/poisson3D-flex.git) and an iterative density calculator (tidal_equilibrium.py).

The first step is to change the "Gitdir" variable in "dir_info.py" to be the local directory where you have cloned the Github repository. This should set all the relevant directories correctly.

Then, you should run "Python LaneEmden.py" to obtain the solution to the LaneEmden equation for the polytropic index "n" that you are interested in. This will produce a spherically symmetric density profile "./data_figs/polytrope_profile.txt", which will be used as the initial guess for density profile iterations under tidal potential.

Then, you should be able to run the code using "Python tidal_equilibrium.py", and if successful, you will see that print-out like the following:

for initial profile, rmax=0.703<br>
qstar(0)=0.30112<br>
Poisson3D iterations: 147<br>
Poisson3D running time: 9.000000<br>
finished iteration 0<br>
qstar(1)=0.30730, frac_delta_q=2.010e-02<br>
Poisson3D iterations: 158<br>
Poisson3D running time: 10.000000<br>
finished iteration 1<br>
qstar(2)=0.30755, frac_delta_q=8.210e-04<br>
Poisson3D iterations: 158<br>
Poisson3D running time: 10.000000<br>
finished iteration 2<br>
qstar(3)=0.31002, frac_delta_q=7.948e-03<br>
....

The subroutines used in the Poisson solver have already been compiled, but if they do not work, you may need to manually compile them. Follow "compile_info.txt" to do that.

It should take <= 10 iterations (and ~10 seconds per iteration for standard resolution Nres=100) for the final density profile to converge and the **final** potential field is saved in "./data_figs/potential%d.txt".

The relevant parameters in the equilibrium tide problem are all in "para_tidal_eq.py", where some recommended values are given. For the case of npoly = 1.5 (corresponding to adiabatic index of gamma=5/3), Kentr=0.5 (entropy constant), Qbh = 1e6 (BH mass in units of Msun), sma=100 (binary separation in units of Rsun), I find that the star marginally fills up its Roche Lobe for a critical central density of rhoc ~ 13 (in units of Msun/Rsun^3).

To visualize the final density/potential profiles, you will run "Python plt_rho_Phi.py a b", which takes two command-line inputs a and b. The first input a is either 0 (plotting the potential field) or 1 (plotting the density field). The second input b specifies the density field at a certain iteration along the way to convergence. For instance, "Python plt_rho_Phi.py 0 7" will plot the final/converged potential field (the L1 point is indicated by a black dot) together with the density field in the 7th iteration. There are also two free parameters in "plt_rho_Phi.py": z_plt specifies the slice at a given z, and Ncont gives the number of contour levels for the potential field. If you run "plt_rho_Phi.py 1 0", "plt_rho_Phi.py 1 1", "plt_rho_Phi.py 1 2", etc, you will see how convergence is reached along the way.
