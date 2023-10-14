# Roche_tidal_equilibrium
The goal is to calculate the density profile of a polytropic star in equilibrium under Newtonian tidal potential of a point-mass companion. The current implementation assumes a star in synchronous rotation and calculations are done in the corotating frame with the coordinate origin at the center of the star.

There are two major pieces of codes: a C-based Poisson solver (poisson3D) from [Derek Harrison's Github repository](https://github.com/derekharrison/poisson3D-flex.git) and an iterative density calculator (tidal_equilibrium.py).

The first step is to change the "dir_main" variable in "dir_info.py" to be the local directory where you have cloned the Github repository. This should set all the relevant directories correctly.

Then, you should run "Python LaneEmden.py" to obtain the solution to the LaneEmden equation for the polytropic index "npoly" (set in para_tidal_eq.py) that you are interested in. This will produce a density profile for a spherically symmetric star at "data_figs/polytrope_profile_npoly%.5f.txt", which will be used as the initial guess for the density profile iterations under tidal potential.

Then, you should be able to run the code with "Python tidal_equilibrium.py 0.578 2.043", where the two numbers are "Kentr" (entropy constant) and "rhoc" (central density) are just an example. This should take a few minutes to finish. The time it takes depends on the resolution, which is set by "Nresz", which is the number of grid points in the z direction, and "Nresx" in the x direction is twice of "Nresz". If successful, you will see output in a log file at data_figs/Kentr0.578/rhoc2.043/output.txt, which says

###################
for initial profile, rmax=2.331
qstar(0)=2.00003
finished iteration 0
warning: no L2 point detected, consider using a larger Lmax!
qstar(0)=2.03977, frac_delta_q=1.948e-02
finished iteration 1
warning: no L2 point detected, consider using a larger Lmax!
qstar(1)=1.99052, frac_delta_q=-2.474e-02
finished iteration 2
warning: no L2 point detected, consider using a larger Lmax!
qstar(2)=2.02434, frac_delta_q=1.671e-02
finished iteration 3
warning: no L2 point detected, consider using a larger Lmax!
qstar(3)=2.01151, frac_delta_q=-6.379e-03
finished iteration 4
warning: no L2 point detected, consider using a larger Lmax!
qstar(4)=2.01432, frac_delta_q=1.397e-03
finished iteration 5
warning: no L2 point detected, consider using a larger Lmax!
qstar(5)=2.01396, frac_delta_q=-1.818e-04
finished iteration 6
warning: no L2 point detected, consider using a larger Lmax!
qstar(6)=2.01397, frac_delta_q=4.173e-06
solution has converged because Mstar(7)=Mstar(6)!
warning: no L2 point detected, consider using a larger Lmax!
equilibrium result: qstar=2.01397, xL1=3.22715, PhitotL1+1.5Q/a=-1.13405320, xsurf=2.82000, Phitotsurf+1.5Q/a=-1.15629677, rhoc=2.04348
FINAL: detached
##################

The warning about "no L2 point detected" is fine, but one could use a large box size which is set by Lmax in "para_tidal_eq.py" to enclose the L2 point. Using a larger box means the boundary conditions are more accurate but the resolution will be poorer.

The subroutines used in the Poisson solver have already been compiled, but if they do not work, you may need to manually compile them. See "poinsson3D/src/compile_info.txt" for more information.

It should take <= 10 iterations (and ~10 seconds per iteration for standard resolution Nresz=100) for the final density profile to converge and the **final** potential field is saved in "data_figs/Kentr%.3f/rhoc%.3f/potential%d.txt". If you set the "OnlySaveLast" parameter to be "True", then only the last=converged results will be saved, meaning that the files generated in intermediate steps will be removed.

The relevant parameters in the equilibrium tide problem are all in "para_tidal_eq.py", where some recommended values are given. For the case of npoly = 1.5 (corresponding to adiabatic index of gamma=5/3), Kentr=0.5 (entropy constant), Qbh = 1e6 (BH mass in units of Msun), sma=100 (binary separation in units of Rsun), I find that the star marginally fills up its Roche Lobe for a critical central density of rhoc ~ 13 (in units of Msun/Rsun^3). 

To visualize the final density/potential profiles, you will run "Python plt_rho_Phi.py", which takes four command-line inputs, as explained in that code. If successful, you should get a figure showing the potential and density profiles of the star! The potential profile is shown by a colormap with black contours, and the density profile is shown by orange contours. The L1 point (and potentially the L2 point) will be shown by an empty black circle, and the stellar surface along the x-axis will be shown by filled black circles.

There are also two more free parameters: z_plt specifies the slice at a given z, and Ncont gives the number of contour levels for the potential field. You can also change the contour levels for the density field by modifying "rhocB_levels".

If you would like to see how the convergence is reached along the iteration steps, you first need to set "OnlySaveLast" to be False. After running tidal_equilibrium.py, you will get all the potential and density profiles for the intermediate steps. You can then run "plt_rho_Phi.py 0 0 Kentr rhoc", "plt_rho_Phi.py 0 1 Kentr rhoc", "plt_rho_Phi.py 0 2 Kentr rhoc", etc to get a series of figures.

Finally, it is often desired to calculate the equilibrium solution for a wide range of entropy constants, and for each entropy constant, we would like to consider a wide range of central densities. For this purpose, you need "run_Kentr_rhoc.py", which does the calculations in parallel. The default case in "run_Kentr_rhoc.py" is for a single Kentr, which is set by a spherically symmetric polytropic star with a given stellar mass (Mstar) and radius (Rstar). For this given Kentr, the code obtains the equilibrium solutions for a wide range of central densities, between "rhocmin" and "rhocmax", in a logarithmic grid with the number of grid points given by "Nrhoc = 2*Ncpu". Here "Ncpu" is the number of processors you want to use.
