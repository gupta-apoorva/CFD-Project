#include "helper.h"
#include "init.h"

int read_parameters( const char *szFileName,   /* name of the file */
                      double *Re,                /* reynolds number   */
                      double *Pr,                /* Prandlt Number*/
                      double *UI,                /* velocity x-direction */
                      double *VI,                /* velocity y-direction */
                      double *PI,                /* pressure */
                      double *TI,
                      double *GX,                /* gravitation x-direction */
                      double *GY,                /* gravitation y-direction */
                      double *t_end,             /* end time */
                      double *xlength,           /* length of the domain x-dir.*/
                      double *ylength,           /* length of the domain y-dir.*/
                      double *dt,                /* time step */
                      double *dx,                /* length of a cell x-dir. */
                      double *dy,                /* length of a cell y-dir. */
                      int  *imax,                /* number of cells x-direction*/
                      int  *jmax,                /* number of cells y-direction*/
                      double *alpha,             /* uppwind differencing factor*/
                      double *omg,               /* relaxation factor */
                      double *tau,               /* safety factor for time step*/
                      double *eps,               /* accuracy bound for pressure*/
                      double *dt_value,           /* time for output */
                      int *wl,
                      int *wr,
                      int *wt,
                      int *wb,
                      double *beta,
                      double *T_l,
                      double *T_r,
                      double *T_t,
                      double *T_b)
{
  
                      READ_DOUBLE( szFileName, *xlength );
                      READ_DOUBLE( szFileName, *ylength );

                      READ_DOUBLE( szFileName, *Re    );
                      READ_DOUBLE( szFileName, *Pr );
                      READ_DOUBLE( szFileName, *t_end );
                      READ_DOUBLE( szFileName, *dt    );

                      READ_INT   ( szFileName, *imax );
                      READ_INT   ( szFileName, *jmax );

                      READ_DOUBLE( szFileName, *omg   );
                      READ_DOUBLE( szFileName, *eps   );
                      READ_DOUBLE( szFileName, *tau   );
                      READ_DOUBLE( szFileName, *alpha );

                      READ_DOUBLE( szFileName, *dt_value );

                      READ_DOUBLE( szFileName, *UI );
                      READ_DOUBLE( szFileName, *VI );
                      READ_DOUBLE( szFileName, *GX );
                      READ_DOUBLE( szFileName, *GY );
                      READ_DOUBLE( szFileName, *PI );
                      READ_DOUBLE( szFileName, *TI );

                      READ_INT( szFileName, *wl);
                      READ_INT( szFileName, *wr);
                      READ_INT( szFileName, *wt);
                      READ_INT( szFileName, *wb);

                      READ_DOUBLE( szFileName, *beta );

                      READ_DOUBLE( szFileName, *T_l );
                      READ_DOUBLE( szFileName, *T_r );
                      READ_DOUBLE( szFileName, *T_t );
                      READ_DOUBLE( szFileName, *T_b );


                      *dx = *xlength / (double)(*imax);
                      *dy = *ylength / (double)(*jmax);

                      return 1;
}

void init_uvpt(double UI,double VI,double PI,double TI, int imax,int jmax,double **U,double **V,double **P, double** T)
{
	init_matrix(U,0,imax+1,0,jmax+1,UI);
	init_matrix(V,0,imax+1,0,jmax+1,VI);
	init_matrix(P,0,imax+1,0,jmax+1,PI);
  init_matrix(T,0,imax+1,0,jmax+1,TI);
}



