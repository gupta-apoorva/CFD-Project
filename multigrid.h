#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include <stdio.h>
#ifndef __PETSCSOLVE_H_
#define __PETSCSOLVE_H_

PetscErrorCode ComputeJacobian(KSP,Mat,Mat,void*);
PetscErrorCode ComputeRHS(KSP,Vec,void*);
PetscErrorCode ComputeTrueSolution(DM, Vec);
PetscErrorCode VecView_VTK(Vec, const char [], const char []);
PetscErrorCode KSPGetSolution(KSP ksp,Vec *v);
typedef enum {DIRICHLET, NEUMANN} BCType;

typedef struct {
  double **RHS;
  BCType      bcType;
} UserContext;

#endif
