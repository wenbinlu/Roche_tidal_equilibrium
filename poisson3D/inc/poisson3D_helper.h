/*
 * poisson3D_helper.h
 *
 *  Created on: Sep 26, 2018
 *      Author: Derek W. Harrison
 */

#ifndef poisson3D_HELPER_H_
#define poisson3D_HELPER_H_

#include "user_types.h"


void vector_addition(double *v1,
                     double f_v1,
                     double *v2,
                     double f_v2,
                     int size_vec,
                     double *vr);

void dot_product(double *v1,
                 double *v2,
                 int size_vec,
                 double *vr);

void mat_vec_mult(grid_size_t grid_size,
                  double** A,
                  double* p,
                  double* Ap);

void generate_coefficient_matrix(domain_size_t domain_size,
                                 grid_size_t grid_size,
                                 boundary_conditions_t boundary_conditions,
                                 physical_paramaters_t physical_parameters,
                                 double (*source)(double x,double y,double z),
                                 grid_coordinates_t* grid_coordinates,
                                 kershaw_algorithm_data_t* kershaw_data);

void Ly_solver(grid_size_t grid_size,
               double** L,
               double* r,
               double* y);

void LTz_solver(grid_size_t grid_size,
                double** L,
                double* y,
                double* z);

void incomplete_cholesky_factorization(grid_size_t grid_size,
                                       double** A,
                                       double** L);

void initialize_temperature_field(grid_size_t grid_size,
                                  double *temp_field);

void generate_grid_coordinates(domain_size_t domain_size,
                               grid_size_t grid_size,
                               grid_coordinates_t* grid_coordinates);

void preconditioning(grid_size_t grid_size,
                     kershaw_algorithm_data_t* kershaw_data);

void execute_kershaw_algorithm(grid_size_t grid_size,
                               kershaw_algorithm_data_t* kershaw_data);

void processing_results(grid_size_t grid_size,
                        kershaw_algorithm_data_t* kershaw_data,
                        double*** T);

void print_results(kershaw_algorithm_data_t* kershaw_data,
                   double time_spent);

#endif /* poisson3D_HELPER_H_ */
