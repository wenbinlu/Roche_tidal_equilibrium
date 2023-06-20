/*
 * boundary_functions.c
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */

#include <stdio.h>
#include <math.h>
#include "../inc/mypara.h"



/*-----------------------------------------------------------------------------------------------*/
double boundary_west(double y, double z)
{
  // NEUMANN
    /*
     * Specify the west face boundary distribution for the Poisson equation
     * equation.
     * If the boundary is of type DIRICHLET, specify the distribution at the boundary. If the boundary
     * is of type NEUMANN specify the flux distribution at the boundary.
     *
     * The distribution or flux distribution is of the form f(y,z).
     *
     * input    y
     * input    z
     *
     * return   f
     */

    double f = 0.0;

    return f;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_east(double y, double z)
{
  // DIRICHLET
    /*
     * Specify the east face boundary distribution for the Poisson equation
     * equation.
     * If the boundary is of type DIRICHLET, specify the distribution at the boundary. If the boundary
     * is of type NEUMANN specify the flux distribution at the boundary.
     *
     * The distribution or flux distribution is of the form f(y,z).
     *
     * input    y
     * input    z
     *
     * return   f
     */
    double x = Lmax;
    double f = -8*Mstar/pow(x*x + y*y + z*z, 0.5);

    return f;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_south(double x, double z)
{
  // NEWMANN
    /*
     * Specify the south face boundary distribution for the Poisson equation
     * equation.
     * If the boundary is of type DIRICHLET, specify the distribution at the boundary. If the boundary
     * is of type NEUMANN specify the flux distribution at the boundary.
     *
     * The distribution or flux distribution is of the form f(x,z).
     *
     * input    x
     * input    z
     *
     * return   f
     */

    double f = 0.0;

    return f;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_north(double x, double z)
{
  // DIRICHLET
    /*
     * Specify the north face boundary distribution for the Poisson equation
     * equation.
     * If the boundary is of type DIRICHLET, specify the distribution at the boundary. If the boundary
     * is of type NEUMANN specify the flux distribution at the boundary.
     *
     * The distribution or flux distribution is of the form f(x,z).
     *
     * input    x
     * input    z
     *
     * return   f
     */

    double y = Lmax;
    double f = -8*Mstar/pow(x*x + y*y + z*z, 0.5);

    return f;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_bottom(double x, double y)
{
    /*
     * Specify the bottom face boundary distribution for the Poisson equation
     * equation.
     * If the boundary is of type DIRICHLET, specify the distribution at the boundary. If the boundary
     * is of type NEUMANN specify the flux distribution at the boundary.
     *
     * The distribution or flux distribution is of the form f(x,y).
     *
     * input    x
     * input    y
     *
     * return   f
     */

    double f = 0.0;

    return f;

}


/*-----------------------------------------------------------------------------------------------*/
double boundary_top(double x, double y)
{
    /*
     * Specify the top face boundary distribution for the Poisson equation
     * equation.
     * If the boundary is of type DIRICHLET, specify the distribution at the boundary. If the boundary
     * is of type NEUMANN specify the flux distribution at the boundary.
     *
     * The distribution or flux distribution is of the form f(x,y).
     *
     * input    x
     * input    y
     *
     * return   f
     */
    double z = Lmax;
    double f = -8*Mstar/pow(x*x + y*y + z*z, 0.5);

    return f;

}



/*-----------------------------------------------------------------------------------------------*/
double source_equation(double x, double y, double z)
{
    /*
     * Specify the source equation for the transient heat conduction equation:
     *
     * gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q = 0
     *
     * The source equation is of the form q(x,y,z).
     *
     * input    x
     * input    y
     * input    z
     *
     * return   q
     */

    //double q = -sin(M_PI*x) * sin(M_PI*y) * sin(M_PI*z);

  double dx = Lmax/Nres;   // grid size
  
  int i = (int) round((x-0.5*dx)/dx);
  int j = (int) round((y-0.5*dx)/dx);
  int k = (int) round((z-0.5*dx)/dx);

  //printf("dx=%f, Mstar=%f", dx, Mstar);
  //printf("rho000=%f", rhoarr[0][0][0]);
  
  return -4*pi * rhoarr[i][j][k];

}
