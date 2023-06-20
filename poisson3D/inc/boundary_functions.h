/*
 * boundary_functions.h
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */

#ifndef BOUNDARY_FUNCTIONS_H_
#define BOUNDARY_FUNCTIONS_H_


double boundary_west(double y, double z);
double boundary_east(double y, double z);
double boundary_south(double x,double z);
double boundary_north(double x,double z);
double boundary_bottom(double x,double y);
double boundary_top(double x,double y);
double source_equation(double x, double y, double z);

#endif /* BOUNDARY_FUNCTIONS_H_ */
