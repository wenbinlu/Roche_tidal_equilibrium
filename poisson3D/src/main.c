/*
 * main.c
 *
 *  Created on: September 25 2018
 *      Author: Derek W. Harrison
 *
 *      This code solves the 3D Poisson equation:
 *
 *      gammax*d2T/dx2 + gammay*d2T/dy2 + gammaz*d2T/dz2 + q(x,y,z) = 0
 *
 *      On a rectangular grid. The code handles mixed boundary conditions, i.e.
 *      both Neumann and Dirichlet boundary conditions
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/mypara.h"   // newly added

#include "../inc/boundary_functions.h"
#include "../inc/export_data.h"
#include "../inc/main_utils.h"
#include "../inc/mappings.h"
#include "../inc/memory_functions.h"
#include "../inc/parameters_and_boundary_types.h"
#include "../inc/poisson3D.h"
#include "../inc/user_types.h"

// global variables
double npoly;  // polytropic index
double pi = M_PI;
double Lmax; // maximum domain size [machine units]
int Nres;  // resolution per dimension
char fdir[] = "../../data_figs/";

double ***rhoarr=NULL;  // density field
double Mstar;  // total stellar mass


/*-----------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    bool                           exportData = FALSE;
    domain_size_t                 domain_size = {0};
    grid_size_t                     grid_size = {0};
    grid_coordinates_t*      grid_coordinates = NULL;
    double***                               T = NULL;
    boundary_conditions_t boundary_conditions = {0};
    boundary_type_faces_t boundary_type_faces = {0};
    physical_paramaters_t physical_parameters = {0};

    // read Lane-Emden profile
    FILE *fp;
    char fname[100];
    int i, j, k;
    char *word;

    int Niter = atoi(argv[1]); // convert the 2nd arg to an integer
    npoly = atof(argv[2]);
    Lmax = atof(argv[3]);
    Nres = atoi(argv[4]);

    double dx = Lmax/Nres;   // grid size
    int Nbuf = (int)(40*Nres);
    char buffer[Nbuf];  // make sure this is long enough for a row

    //printf("%i, %f, %f, %i\n", Niter, npoly, Lmax, Nres);

    rhoarr = matrix3D(Nres, Nres, Nres); // allocate memory

    sprintf(fname, "%s%s%i%s", fdir, "rho", Niter, ".txt");
    //printf("%s\n", fname);
    
    fp = fopen(fname, "r");
    
    for (i=0; i<Nres; i++) {
      for (j=0; j<Nres; j++) {
	fgets(buffer, sizeof(buffer), fp);
	//if (i+j == 0) { printf("%s", buffer); }
	k = 0;
	word = strtok(buffer, "\t");
	while (word != NULL) {
	  rhoarr[i][j][k] = atof(word);
	  word = strtok(NULL, "\t");
	  k++;
	}	
      }
    }

    // calculate total stellar mass
    Mstar = 0.0;
    for (i=0; i<Nres; i++) {
      for (j=0; j<Nres; j++) {
	for (k=0; k<Nres; k++) {
	  Mstar = Mstar + rhoarr[i][j][k] * pow(dx, 3);
	}
      }
    }
    //printf("stellar mass = %f\n", Mstar);
    
    
    /* Set parameters and boundary conditions */
    set_parameters_and_boundary_types(&domain_size,
                                      &grid_size,
                                      &physical_parameters,
                                      &boundary_type_faces,
                                      &exportData);
    
    /* Allocating memory for input and output of poisson solver */
    grid_coordinates = allocate_mem_grid_coordinates(grid_size.nx+1, grid_size.ny+1, grid_size.nz+1);
    T                = matrix3D(grid_size.nx+1, grid_size.ny+1, grid_size.nz+1);

    
    /* Setting boundary types and conditions (after reading LaneEmden profile) */
    set_boundary_conditions(boundary_type_faces,
                            &boundary_conditions);

    
    /* Calling Poisson solver */
    poisson3D(domain_size,
              grid_size,
              boundary_conditions,
              physical_parameters,
              &source_equation,
              grid_coordinates,
              T);

    

    /* Exporting data */
    if(exportData)
    {
        export_data(grid_size, grid_coordinates, T);
    }


    /* Freeing memory */
    free_grid_coordinates(grid_coordinates, grid_size.nx+1, grid_size.ny+1);
    free_memory_3D(T, grid_size.nx+1, grid_size.ny+1);
    
    free_memory_3D(rhoarr, Nres, Nres); // added
    
    return 0;

}

