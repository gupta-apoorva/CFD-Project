#include "boundary_val.h"
#include <string.h>
#include <stdio.h>
#include "math.h"

#define NO_SLIP 1         // defining the different kinds of boundary  conditions
#define FREE_SLIP 2
#define OUTFLOW 3


/*#define B_N 2          
#define B_S 1             // Defining the values to find boundary directions. 
#define B_W 4
#define B_E 8
#define B_NE 9
#define B_NW 5
#define B_SE 10
#define B_SW 6*/

void boundaryvalues(int imax,int jmax,int wl, int wr, int wt, int wb, double TI, double T_body,  double **U,double **V, 
	                double** P, double** T, double** G, double** F, int** FLAG, double T_l, double T_r, double T_t, double T_b)
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
	for(int i=1; i<=imax; i++){
        	P[i][0] = P[i][1];
			P[i][jmax+1] = P[i][jmax];
			T[i][0] = 2*T_b - T[i][1];
			T[i][jmax+1] = 2*T_t - T[i][jmax];
        	G[i][0] = V[i][0];
        	G[i][jmax] = V[i][jmax]; 
    	}
    
	for(int j=1; j<=jmax; j++){
        	P[0][j] = P[1][j];
			P[imax+1][j] = P[imax][j];
        	T[0][j] = 2*T_l + T[1][j];
			T[imax+1][j] = 2*T_r + T[imax][j];
        	F[0][j] = U[0][j];
        	F[imax][j] = U[imax][j];
	}


// Taking care of the boundary generated by the arbitrary boundary and all of them are no slip condition...
	/*for (int i = 1; i <=imax; ++i)
	{	
		for (int j = 1; j <=jmax; ++j)
		{
			if (FLAG[i][j] >=1 && FLAG[i][j] <=15)
			{
				if (FLAG[i][j] == B_E)
				{
					U[i][j] = 0;
					V[i][i-1] = -V[i+1][j-1];
					V[i][j] = -V[i+1][j];
					F[i][j] = U[i][j];
					P[i][j] = P[i+1][j];
					T[i][j] = 2*T_body - T[i+1][j];
				}
				else if (FLAG[i][j] == B_W)
				{
					U[i-1][j] = 0;
					V[i][j-1] = -V[i-1][j-1];
					V[i][j] = -V[i-1][j];
					F[i-1][j] = U[i-1][j];
					P[i][j] = P[i-1][j];
					T[i][j] = 2*T_body - T[i-1][j];
				}
				else if (FLAG[i][j]  == B_S)
				{
					V[i][j-1] = 0;
					U[i][j] = -U[i][j-1];
					U[i-1][j] = -U[i-1][j-1];
					G[i][j-1] = V[i][j-1];
					P[i][j] = P[i][j-1];
					T[i][j] = 2*T_body - T[i][j-1]; 
				}
				else if (FLAG[i][j]  == B_N)
				{
					V[i][j] =0;
					U[i-1][j] = -U[i-1][j+1];
					U[i][j] = -U[i][j+1];
					G[i][j] = V[i][j];
					P[i][j] = P[i][j+1];
					T[i][j] = 2*T_body - T[i][j+1];				
				}
				else if (FLAG[i][j] == B_NE)
				{
					U[i][j] = 0;
					V[i][j] = 0;
					U[i-1][j] = -U[i-1][j+1];
					V[i][j-1] = -V[i+1][j-1];
					F[i][j] = U[i][j];
					G[i][j] = V[i][j];
					P[i][j] = (P[i][j+1]+P[i+1][j])/2;
					T[i][j] = 2*T_body - (T[i][j+1]+T[i+1][j])/2;
				}
				else if (FLAG[i][j] == B_NW)
				{
					U[i-1][j] = 0;
					V[i][j] = 0;
					U[i][j] = -U[i][j+1];
					V[i][j-1]= -V[i-1][j-1];
                    F[i-1][j] = U[i-1][j];
                   	G[i][j] = V[i][j];
                   	P[i][j] = (P[i][j+1] + P[i-1][j])/2;
                   	T[i][j] = 2*T_body - (T[i][j+1] + T[i-1][j])/2;
				}
				else if (FLAG[i][j] == B_SE)
				{
					U[i][j] = 0;
					V[i][j-1] = 0;
					U[i-1][j] = - U[i-1][j-1];
					V[i][j] = -V[i+1][j];
					F[i][j] = U[i][j];
					G[i][j-1] = V[i][j-1];
					P[i][j] = (P[i+1][j] + P[i][j-1])/2;
					T[i][j] = 2*T_body - (T[i+1][j] + T[i][j-1])/2;
				}
				else if (FLAG[i][j] == B_SW)
				{
					U[i-1][j] = 0;
					V[i][j-1] = 0;
					U[i][j] = - U[i][j-1];
					V[i][j] = -V[i-1][j];
					F[i-1][j] = U[i][j-1];
					G[i][j-1] = V[i][j-1];
					P[i][j] = (P[i-1][j] + P[i][j-1])/2;	
					T[i][j] = 2*T_body - (T[i-1][j] + T[i][j-1])/2;			
				}

			}
		}
	}*/

}
