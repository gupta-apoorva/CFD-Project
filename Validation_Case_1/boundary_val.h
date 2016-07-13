#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int imax,int jmax,int wl, int wr, int wt, int wb, double TI, double **U,double **V, double** P, 
					double** T, double** G, double** F, double T_l, double T_r, double T_t, double T_b);


#endif
