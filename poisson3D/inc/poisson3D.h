/*
 * poisson3D.h
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#ifndef poisson3D_H_
#define poisson3D_H_

#include "user_types.h"


void poisson3D(domain_size_t domain_size,
               grid_size_t grid_size,
               boundary_conditions_t boundary_conditions,
               physical_paramaters_t physical_parameters,
               double (*source)(double x,double y,double z),
               grid_coordinates_t* grid_coordinates,
               double ***T);

#endif /* poisson3D_H_ */
