#include "boundary_val.h"
#include <string.h>
#include <stdio.h>
#include "math.h"
#include <omp.h>

#define NO_SLIP 1         // defining the different kinds of boundary  conditions
#define FREE_SLIP 2
#define OUTFLOW 3

void boundaryvalues(int imax,int jmax,int wl, int wr, int wt, int wb, double TI, double **U,double **V, double** P, 
					double** T, double** G, double** F, double T_l, double T_r, double T_t, double T_b)

{
// Taking care of the left boundary of the domain	

	if (wl==NO_SLIP)
	{
		for (int j = 1; j <=jmax ; ++j)                           
		{
			V[0][j] = -V[1][j];
			U[0][j] = 0;
		}
	}
	else if (wl==FREE_SLIP)
	{
		for (int j = 1; j <=jmax ; ++j)
		{
			V[0][j] = V[1][j];
			U[0][j] = 0;
		}
	}
	else if (wl==OUTFLOW)
	{
		for (int j = 1; j <=jmax ; ++j)
		{
			V[0][j] = V[1][j];
			U[0][j] = U[1][j];
		}
	}


// Taking care of the right boundary of the domain
	if (wr==NO_SLIP){
		for (int j = 1; j <=jmax ; ++j){
			V[imax+1][j] = -V[imax][j];
			U[imax][j] = 0;
		}
	}
	else if (wr==FREE_SLIP){
		for (int j = 1; j <=jmax ; ++j){
			V[imax+1][j] = V[imax][j];
			U[imax][j] = 0;
		}
	}
	else if (wr==OUTFLOW){
		for (int j = 1; j <=jmax ; ++j){
			V[imax+1][j] = V[imax][j];
			U[imax][j] = U[imax-1][j];
		}
	}


// Taking care of the top boundary of the domain
	if (wt==NO_SLIP){
		for (int i = 1; i <=imax ; ++i){
			U[i][jmax+1] = -U[i][jmax];
			V[i][jmax] = 0;
		}
	}
	else if (wt==FREE_SLIP){
		for (int i = 1; i <=imax ; ++i)	{
			U[i][jmax+1] = U[i][jmax];
			V[i][jmax] = 0;
		}
	}
	else if (wt==OUTFLOW){
		for (int i = 1; i <=imax ; ++i){
			U[i][jmax+1] = U[i][jmax];
			V[i][jmax] = V[i][jmax-1];
		}
	}


// Taking care of the bottom boundary of the domain
	if (wb==NO_SLIP){
		for (int i = 1; i <=imax ; ++i){
			U[i][0] = -U[i][1];
			V[i][0] = 0;
	        }
        }

	else if (wb==FREE_SLIP){
		for (int i = 1; i <=imax ; ++i)	{
			U[i][0] = U[i][1];
			V[i][0] = 0;
	         }
        }

	else if (wb==OUTFLOW){
		for (int i = 1; i <=imax ; ++i){
			U[i][0] = U[i][1];
			V[i][0] = V[i][1];
	         }
       }

// Since the pressure P , T , F and G dosen't depend on the type of boundary so taking care of them together....
    #pragma omp parallel for
	for(int i=1; i<=imax; i++){
        	P[i][0] = P[i][1];
			P[i][jmax+1] = P[i][jmax];
			T[i][0] = 2*T_b - T[i][1];
			T[i][jmax+1] = 2*T_t - T[i][jmax];
        	G[i][0] = V[i][0];
        	G[i][jmax] = V[i][jmax]; 
    	}
    #pragma omp parallel for
	for(int j=1; j<=jmax; j++){
        	P[0][j] = P[1][j];
			P[imax+1][j] = P[imax][j];
        	T[0][j] = 2*T_l + T[1][j];
			T[imax+1][j] = 2*T_r + T[imax][j];
        	F[0][j] = U[0][j];
        	F[imax][j] = U[imax][j];
	}

}
