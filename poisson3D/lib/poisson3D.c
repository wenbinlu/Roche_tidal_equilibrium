/*
 * poisson3D.c
 *
 *  Created on: Sep 25, 2018
 *      Author: Derek W. Harrison
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../inc/poisson3D_helper.h"
#include "../inc/memory_functions.h"
#include "../inc/user_types.h"


/*-----------------------------------------------------------------------------------------------*/
void poisson3D(domain_size_t domain_size,
               grid_size_t grid_size,
               boundary_conditions_t boundary_conditions,
               physical_paramaters_t physical_parameters,
               double (*source)(double x,double y,double z),
               grid_coordinates_t* grid_coordinates,
               double ***T)
{
    /*
     * This function solves the 3D modified poisson equation:
     *
     * gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q(x,y,z) = 0
     *
     * on a rectangluar grid. The code handles mixed boundary conditions, i.e.
     * both Neumann and Dirichlet boundary conditions.
     *
     * The equation is solved using a preconditioned conjugate gradient method with
     * incomplete cholesky factorization as preconditioner.
     *
     * Note: a node numbering scheme was used such that the first element of all
     * arrays start at 1 as opposed to zero. The numbering scheme used in this work is
     *
     * nn = i + (j-1)*nx + (k-1)*nx*ny
     *
     * where nn is the node number, nx the number of nodes used in the x direction,
     * ny the number of nodes used in the y direction and i, j, k are indices ranging
     * from 1 to nt where nt is the total number of nodes in the 3D grid, i.e.: nt = nx*ny*nz
     * with nz the number of nodes in the z direction.
     *
     * input     domain_size
     * input     grid_size
     * input     boundary_conditions
     * input     physical_parameters
     * input     source()
     * output    grid_coordinates
     * output    T
     */

    clock_t begin, end;
    double time_spent;
    kershaw_algorithm_data_t* kershaw_data = NULL;


    /* Allocating memory */
    kershaw_data = allocate_kershaw_data(grid_size);


    /* Initialize data */
    initialize_temperature_field(grid_size,
                                 kershaw_data->x);

    generate_grid_coordinates(domain_size,
                              grid_size,
                              grid_coordinates);

    /* Execute kershaw algorithm */
    begin = clock();

    generate_coefficient_matrix(domain_size,
                                grid_size,
                                boundary_conditions,
                                physical_parameters,
                                source,
                                grid_coordinates,
                                kershaw_data);

    preconditioning(grid_size,
                    kershaw_data);

    execute_kershaw_algorithm(grid_size,
                              kershaw_data);

    end = clock();


    /* Process results */
    processing_results(grid_size,
                       kershaw_data,
                       T);


    /* Print some results */
    time_spent = (end - begin)/CLOCKS_PER_SEC;
    //print_results(kershaw_data, time_spent);


    /* Deallocate data */
    free_kershaw_data(kershaw_data, grid_size);

}
