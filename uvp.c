#include "uvp.h"
#include "helper.h"
#include <stdlib.h>
#include "math.h"

/* calculate the tine stepping based on the value of tau*/

void calculate_dt(
  double Re,
  double Pr,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V)

{
if (tau>0)
  {
    double umax = 0;
    double vmax = 0;
    for (int i=1;i<=imax;i++) 
    {
      for (int j=1;j<=jmax;j++)
      {
         if (abs(U[i][j]) > umax)
            {
               umax = abs(U[i][j]);
            }
        
       }
    }
    for (int i=1;i<=imax;i++) 
    {
      for (int j=1;j<=jmax;j++)
      {
          if (abs(V[i][j]) > vmax)
            {
               vmax = abs(V[i][j]);
            }      
       }
    }
    
    *dt = tau*fmin(fmin(dx/umax,dy/vmax),fmin(Re/(2*(1/(dx*dx)+1/(dy*dy))),(Re*Pr)/(2*(1/(dx*dx)+1/(dy*dy)))));
   }

}


/* Calculating right hand size of the pressure equation */

void calculate_rs(double dt,double dx,double dy,  int imax,  int jmax,  double **F,  double **G,  double **RS,int** FLAG)
{
for (int i=1;i<=imax;i++)
   {
     for (int j=1;j<=jmax;j++)
      {
        if (FLAG[i][j] >=16 && FLAG[i][j] <=31)
         RS[i][j] = 1/dt*((F[i][j] - F[i-1][j])/dx + (G[i][j] - G[i][j-1])/dy);
      }
   }
}


/* updating the values of U and V */

void calculate_uv(double dt,  double dx,  double dy,  int imax,  int jmax,double **U,  double **V,  double **F, double **G,double **P,int** FLAG)
{
for (int i=1;i<=imax-1;i++)
  {
    for (int j =1;j<=jmax;j++)
      {
        if (FLAG[i][j] >=16 && FLAG[i][j] <=31)
        U[i][j] = F[i][j] - dt/dx*(P[i+1][j] - P[i][j]);
      }
   }

for (int i=1;i<=imax;i++)
  {
    for (int j =1;j<=jmax-1;j++)
      {
        if (FLAG[i][j] >=16 && FLAG[i][j] <=31)
        V[i][j] = G[i][j] - dt/dy*(P[i][j+1] - P[i][j]);
      }
  }

}

// calculate_fg function
 
void calculate_fg(double Re,double GX,double GY,double alpha,double dt,double dx,  double dy,int imax,int jmax, double beta, double **U,  double **V, double**T, double **F,  double **G, int** FLAG)
  {
   double duvdy,du2dx,d2udx2,d2udy2,d2vdx2,d2vdy2,duvdx,dv2dy;
   for (int i = 1; i <= imax-1; ++i)
  {
    for (int j = 1; j <= jmax; ++j)
    {
      if (FLAG[i][j] >=16 && FLAG[i][j] <=31)
      {
        duvdy = 1/(4*dy)*(((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1]))-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])+alpha*((abs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1]))-abs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])));
        du2dx = 1/(4*dx)*(((pow((U[i][j]+U[i+1][j]),2))-(pow((U[i-1][j]+U[i][j]),2)))+alpha*((abs(U[i][j]+U[i+1][j]))*(U[i][j]-U[i+1][j])-(abs(U[i-1][j]+U[i][j]))*(U[i-1][j]-U[i][j])));
        d2udx2 = (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
        d2udy2 = (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
       
        F[i][j] = U[i][j] + dt*(1/Re*(d2udx2 + d2udy2) - du2dx - duvdy + GX) - beta*dt/2*(T[i][j] + T[i+1][j])*GX;
      }
    }
  }
  for (int i =1;i<=imax;i++)
  {
    for (int j =1;j<=jmax-1;j++)
    {
      if (FLAG[i][j] >=16 && FLAG[i][j] <=31)
      {
        d2vdx2 = (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
        d2vdy2 = (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
        duvdx = 1/(4*dx)*(((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j]))-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])+alpha*((abs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j]))-abs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])));
        dv2dy = 1/(4*dy)*(pow(((V[i][j] + V[i][j+1])),2) - pow(((V[i][j-1] + V[i][j])),2) + alpha*(abs(V[i][j] + V[i][j+1])*(V[i][j] - V[i][j+1]) - abs(V[i][j-1] + V[i][j])*(V[i][j-1] - V[i][j])));
      
        G[i][j] = V[i][j] + dt*(1/Re*(d2vdx2 + d2vdy2) - duvdx - dv2dy + GY) - beta*dt/2*(T[i][j] + T[i][j+1])*GY;
      }
    }
   }
  }

void update_T(int imax, int jmax, double Re, double Pr, double dt, double dx, double dy, double alpha, double** T, double** U, double** V, int** FLAG)
  {
  double dUTdx, dVTdy, d2Tdx2, d2Tdy2;
  for (int i =1;i<=imax;i++)
  {
   for (int j =1;j<=jmax-1;j++)
    {
      if (FLAG[i][j] >=16 && FLAG[i][j] <=31)
      {
        dUTdx = (U[i][j]*(T[i][j] + T[i+1][j])/2 - U[i-1][j]*(T[i-1][j] + T[i][j])/2)/dx + (abs(U[i][j])*(T[i][j] - T[i+1][j])/2 - abs(U[i-1][j])*(T[i-1][j] - T[i][j])/2)*alpha/dx;
        dVTdy = (V[i][j]*(T[i][j] + T[i][j+1])/2 - V[i][j-1]*(T[i][j-1] + T[i][j])/2)/dy + (abs(V[i][j])*(T[i][j] - T[i][j+1])/2 - abs(V[i][j-1])*(T[i][j-1] - T[i][j])/2)*alpha/dy;
        d2Tdx2 = (T[i+1][j] - 2*T[i][j] + T[i-1][j])/pow(dx,2);
        d2Tdy2 = (T[i][j+1] - 2*T[i][j] + T[i][j-1])/pow(dy,2);

        T[i][j] = T[i][j] + dt*(1/(Re*Pr)*(d2Tdx2 + d2Tdy2) - dUTdx - dVTdy);
      }
    }
  }
  }
