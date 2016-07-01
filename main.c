#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>
#include <string.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_FLAG(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

#define TILTED_PLATE 1    // Defining different problem types
#define PLANE_SHEAR 2
#define FLOW_STEP 3
#define C_F 1
#define C_B 0

int main(int argn, char** args)
{

   double** U;
   double** V;
   double** P;
   double** T;
   double Re;
   double Pr;
   double beta;               
   double UI;                
   double VI;               
   double PI; 
   double TI;          
   double GX;                
   double GY;                
   double t_end;            
   double xlength;          
   double ylength;           
   double dt;                
   double dx;             
   double dy;               
   int  imax;                
   int  jmax;               
   double alpha;            
   double omg;               
   double tau;              
   int  itermax;             
   double eps;              
   double dt_value;          
   double** RS;
   double** F;
   double** G;
   int** FLAG;
   int** pgm = NULL;
   int wl;	
   int wr;
   int wt;
   int wb;
   double T_body;
   double T_l;
   double T_r;
   double T_t;
   double T_b;
   
	
//setting the parameters
	read_parameters( "problem.dat", &Re ,&Pr, &UI , &VI, &PI, &TI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax,
		         &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &wl, &wr, &wt, &wb, &beta, &T_body,&T_l, &T_r, &T_t, &T_b);

	  pgm = read_pgm("mesh.pgm");
   printf("hello " );

// Creating the arrays U,V and P
	  U = matrix ( 0 , imax+1 , 0 , jmax+1 );
	  V = matrix ( 0 , imax+1 , 0 , jmax+1 );
	  P = matrix ( 0 , imax+1 , 0 , jmax+1 );
    T = matrix ( 0 , imax+1 , 0 , jmax+1 );
        

// Creating arrays for right side of pressure poissons equation (RS) and F and G
	  RS = matrix ( 0, imax+1, 0, jmax+1);
	  F =  matrix ( 0, imax+1, 0, jmax+1);
	  G =  matrix ( 0, imax+1, 0, jmax+1);
	  FLAG = imatrix ( 0, imax+1, 0, jmax+1);

// Initializing the arrays U,V,P,RS,F, G and FLAG field
	  init_uvpt( UI, VI,PI,TI,imax, jmax,U,V,P,T);
	  init_matrix(RS, 0, imax+1, 0, jmax+1, 0);
	  init_matrix(F, 0, imax+1, 0, jmax+1, 0);
  	init_matrix(G, 0, imax+1, 0, jmax+1, 0);
  	init_imatrix(FLAG, 0, imax+1, 0, jmax+1, C_F);

// Setting the value to be C_b by checking whether it is a fluid or boundary...
	  for (int i = 0; i < imax+2; ++i){
	  for (int j = 0; j < jmax+2; ++j){
	   	if (pgm[i][j] == 0)  
  		FLAG[i][j] = C_B;	
  	}
  	}

// Making a dummy matrix as it will be required later to set the flag field...

  	int dummy[imax+4][jmax+4];

  	for (int j = 0; j < jmax+4; ++j){
  	for (int i = 0; i < imax+4; ++i)
        dummy[i][j] = 0;
        }

  	for (int j = 0; j < jmax+2; ++j){
  	for (int i = 0; i < imax+2; ++i)
        dummy[i+1][j+1] = FLAG[i][j];
        }

// Setting the proper values for our flagfield based on the boundary..

  	for (int j = 0; j < jmax+2; ++j){                
  	for (int i = 0; i < imax+2; ++i)
        FLAG[i][j] = 16*dummy[i+1][j+1] + 8*dummy[i+2][j+1] + 4*dummy[i][j+1] + 2*dummy[i+1][j] + dummy[i+1][j+2]; 				
  	}

   /* for (int i = 0 ; i<imax+2 ; i++)
    {
      for (int j = 0 ; j< jmax+1 ; j++)
      {
        printf ("%d ", FLAG[i][j]);
      }
      printf("\n");
    }*/

    double t=0;   // initialize the time
  	int n = 0;    // number of time steps



int count = 0;
while (t<t_end)
  {

// Calculating the proper time step to maintain stability...

      calculate_dt(Re,Pr,tau,&dt,dx,dy,imax,jmax,U,V);    

// Setting the boundary values for U,V,P depending on the flagfield...

      boundaryvalues(imax, jmax, wl , wr, wt, wb , TI, T_body, U, V , P, T, G, F, FLAG, T_l, T_r, T_t, T_b);

    /*      for (int i = 0 ; i<imax+2 ; i++)
    {
      for (int j = 0 ; j< jmax+1 ; j++)
      {
        printf ("%f ", T[i][j]);
      }
      printf("\n");
    }*/

// Setting the special boundary values based on the problem type...

      //spec_boundary_val (pType, imax, jmax, U, V,delta_p,input_vel, Re, ylength, xlength);

// Updating the T matrix for the nxt time step.
      update_T(imax, jmax, Re, Pr, dt, dx, dy, alpha, T, U, V, FLAG);

// Calculating the F and G matrix for only the fluid cells based on the flagfield...

      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, imax, jmax, beta, U, V, T, F, G, FLAG);

// calculating the right hand side of the pressure poissons for only fluid cells based on the flagfield...

      calculate_rs(dt,dx,dy, imax,jmax, F, G, RS,FLAG);

// iteration number...

      int it = 0;

// residual... 
  
      double res = 1000;  

// Doing the successive over relaxation...
      while(it<itermax && res > eps) 
          {
            sor(omg, dx,dy,imax,jmax, P, RS, &res,FLAG);
            it++; 
          }

// calculating the U and V velocities matrices...

      calculate_uv(dt,dx, dy,imax,jmax,U,V,F,G,P,FLAG);

//Going to the next time step...

      t = t+dt;
      n = n+1;  

      if (count % (int)dt_value == 0)
      write_vtkFile("szProblem.vtk", n, xlength, ylength, imax, jmax,dx, dy, U, V, P, T, FLAG);
      
      count ++;
  }

  //write_vtkFile("szProblem.vtk", n, xlength, ylength, imax, jmax,dx, dy, U, V, P, T, FLAG);

free_matrix(U,0,imax+1,0,jmax+1);
free_matrix(V,0,imax+1,0,jmax+1);
free_matrix(P,0,imax+1,0,jmax+1);
free_matrix(T,0,imax+1,0,jmax+1);
free_matrix(RS,0,imax+1,0,jmax+1);
free_matrix(F,0,imax+1,0,jmax+1);
free_matrix(G,0,imax+1,0,jmax+1);
free_imatrix(FLAG,0,imax+1,0,jmax+1);
free_imatrix(pgm,0,imax,0,jmax);

return 0;
}
