/*
 * parameters_and_boundary_types.c
 *
 *  Created on: Oct 25, 2018
 *      Author: derek
 */

#include <math.h>
#include <stdio.h>
#include "../inc/user_types.h"
#include "../inc/mypara.h"


/*-----------------------------------------------------------------------------------------------*/
void set_parameters_and_boundary_types(domain_size_t* domain_size,
                                       grid_size_t* grid_size,
                                       physical_paramaters_t* physical_parameters,
                                       boundary_type_faces_t* boundary_type_faces,
                                       bool* exportData)
{
    /*
     * Specify parameters and boundary types for 3D poisson solver
     *
     * output   domain_size
     * output   grid_size
     * output   physical_parameters
     * output   boundary_type_faces
     * output   exportData
     */
  
  /* Set parameters and boundary conditions */
    domain_size->Lx = 2*Lmax;                               //length of domain along x coordinate
    domain_size->Ly = Lmax;                               //length of domain along y coordinate
    domain_size->Lz = Lmax;                               //length of domain along z coordinate

    // setting the resolution
    grid_size->nx = 2*Nresz;                                 //amount of nodes along x coordinate
    grid_size->ny = Nresz;                                   //amount of nodes along y coordinate
    grid_size->nz = Nresz;                                   //amount of nodes along z coordinate

    // fixing these three
    physical_parameters->conductivity.gammax = 1.0;     //conductivity along x coordinate
    physical_parameters->conductivity.gammay = 1.0;     //conductivity along y coordinate
    physical_parameters->conductivity.gammaz = 1.0;     //conductivity along z coordinate
    
    // (DIRICHLET vs. NEUMANN)  
    boundary_type_faces->west_boundary   = DIRICHLET;    //west face boundary type (yz plane)
    boundary_type_faces->east_boundary   = DIRICHLET;    //east face boundary type
    boundary_type_faces->south_boundary  = NEUMANN;    //south face boundary type (xz plane)
    boundary_type_faces->north_boundary  = DIRICHLET;    //north face boundary type
    boundary_type_faces->bottom_boundary = NEUMANN;    //bottom face boundary type (xy plane)
    boundary_type_faces->top_boundary    = DIRICHLET;    //top face boundary type

    *exportData = TRUE;                                  //export data guard

}
