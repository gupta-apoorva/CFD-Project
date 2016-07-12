///////// HOW TO RUN /////////////////

    //./sim -da_grid_x 65 -da_grid_y 17 -pc_type mg -pc_mg_levels 1 -mg_levels_0_pc_type lu -mg_levels_0_pc_factor_shift_type NONZERO -ksp_monitor
    //./sim -da_grid_x 65 -da_grid_y 17 -pc_type mg -da_refine 0 -ksp_atol 0.0001 -ksp_monitor
    //./sim -da_grid_x 33 -da_grid_y 9 -pc_type mg -da_refine 1 -ksp_atol 0.001 -ksp_monitor
    //./sim -da_grid_x 17 -da_grid_y 5 -pc_type mg -da_refine 2 -ksp_atol 0.001 -ksp_monitor
    //./sim -da_grid_x 9 -da_grid_y 3 -pc_type mg -da_refine 3 -ksp_atol 0.001 -ksp_monitor



#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include <stdio.h>
#include <string.h>
#include "multigrid.h"

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
    int** pgm = NULL;
    int wl;	
    int wr;
    int wt;
    int wb;
    double T_l;
    double T_r;
    double T_t;
    double T_b;
    KSP            ksp;
    DM             da;
    UserContext    user;
    PetscInt       bc;
    PetscErrorCode ierr;
    PetscScalar    **new;
    Vec		 x;
	
//setting the parameters
  	read_parameters( "problem.dat", &Re ,&Pr, &UI , &VI, &PI, &TI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax,
  		         &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &wl, &wr, &wt, &wb, &beta, &T_l, &T_r, &T_t, &T_b);

  	pgm = read_pgm("mesh.pgm");

// Initializing PETSC by giving arguments in the command line.
    PetscInitialize(&argn,&args,(char*)0,NULL);
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,-11,-11,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da);CHKERRQ(ierr);
    ierr = KSPSetDM(ksp,(DM)da);
    ierr = DMSetApplicationContext(da,&user);CHKERRQ(ierr);
    bc          = (PetscInt)NEUMANN; 
    user.bcType = (BCType)bc;

// Calculating the Jocobian needed to solve pressure poissons
    ierr = KSPSetComputeOperators(ksp,ComputeJacobian,&user);CHKERRQ(ierr);

// It basically set the option for the solver from the command line.    
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

// Creating the arrays U,V and P
	  U = matrix ( 0 , imax+1 , 0 , jmax+1 );
	  V = matrix ( 0 , imax+1 , 0 , jmax+1 );
	  P = matrix ( 0 , imax+1 , 0 , jmax+1 );
    T = matrix ( 0 , imax+1 , 0 , jmax+1 );
        

// Creating arrays for right side of pressure poissons equation (RS) and F and G
	  RS = matrix ( 0, imax+1, 0, jmax+1);
	  F =  matrix ( 0, imax+1, 0, jmax+1);
	  G =  matrix ( 0, imax+1, 0, jmax+1);

// Initializing the arrays U,V,P,RS,F, G field
	  init_uvpt( UI, VI,PI,TI,imax, jmax,U,V,P,T);
	  init_matrix(RS, 0, imax+1, 0, jmax+1, 0);
	  init_matrix(F, 0, imax+1, 0, jmax+1, 0);
  	init_matrix(G, 0, imax+1, 0, jmax+1, 0);

// initialize the time
    double t=0; 

// number of time steps      
  	int n = 0;    

    int count = 0; 
    while (t<t_end)
    {

// Calculating the proper time step to maintain stability...

      calculate_dt(Re,Pr,tau,&dt,dx,dy,imax,jmax,U,V);    

// Setting the boundary values for U,V,P

      boundaryvalues(imax, jmax, wl , wr, wt, wb , TI, U, V , P, T, G, F, T_l, T_r, T_t, T_b);

// Updating the T matrix for the nxt time step.
      update_T(imax, jmax, Re, Pr, dt, dx, dy, alpha, T, U, V);

// Calculating the F and G matrix 

      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, imax, jmax, beta, U, V, T, F, G);

// calculating the right hand side of the pressure poissons

      calculate_rs(dt,dx,dy, imax,jmax, F, G, RS);

      user.RHS  = RS;

// Computing the RHS as PETSC wants
      ierr = KSPSetComputeRHS(ksp,ComputeRHS,&user);CHKERRQ(ierr);

// Solving the prssure poissons equation
      ierr = KSPSolve(ksp,NULL,NULL);CHKERRQ(ierr);

// Extracting the solution form the solver
      KSPGetSolution(ksp,&x);
      ierr =DMDAVecGetArray(da, x, &new);CHKERRQ(ierr);

// Putting the values from the solver at correct location in P array.
      for(int i = 0;i<imax;i++){
      for(int j = 0;j<jmax;j++){
         P[i+1][j+1] =new[j][i];
      }
      }

// calculating the U and V velocities matrices...

      calculate_uv(dt,dx, dy,imax,jmax,U,V,F,G,P);

//Going to the next time step...

      t = t+dt;
      n = n+1;

      printf("time = %f",t);  

      if (count % (int)dt_value == 0)
      write_vtkFile("szProblem.vtk", n, xlength, ylength, imax, jmax,dx, dy, U, V, P, T);
      
      count ++;
    }

  //write_vtkFile("szProblem.vtk", n, xlength, ylength, imax, jmax,dx, dy, U, V, P, T);

free_matrix(U,0,imax+1,0,jmax+1);
free_matrix(V,0,imax+1,0,jmax+1);
free_matrix(P,0,imax+1,0,jmax+1);
free_matrix(T,0,imax+1,0,jmax+1);
free_matrix(RS,0,imax+1,0,jmax+1);
free_matrix(F,0,imax+1,0,jmax+1);
free_matrix(G,0,imax+1,0,jmax+1);
free_imatrix(pgm,0,imax,0,jmax);
ierr = DMDestroy(&da);CHKERRQ(ierr);
ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}
