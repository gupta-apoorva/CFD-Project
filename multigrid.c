#include "multigrid.h"
 
PetscErrorCode ComputeRHS(KSP ksp,Vec b,void *ctx)
{
  UserContext    *user = (UserContext*)ctx;
  PetscErrorCode ierr;
  PetscInt       i,j,M,N,xm,ym,xs,ys;
  PetscScalar    **array;
  DM             da;
  PetscScalar    Hx,Hy;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da, 0, &M, &N, 0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  Hx   = 1.0/(PetscReal)(M);
  Hy   = 1.0/(PetscReal)(N);

  ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, b, &array);CHKERRQ(ierr);


  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      array[j][i] =-Hx*Hy*user->RHS[i+1][j+1];
    }
  }

  ierr = DMDAVecRestoreArray(da, b, &array);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);


 /* if (user->bcType == NEUMANN) {
    MatNullSpace nullspace;

    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
    ierr = MatNullSpaceRemove(nullspace,b);CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
  }*/
  PetscFunctionReturn(0);
}


PetscErrorCode ComputeJacobian(KSP ksp,Mat J, Mat jac,void *ctx)
{
  UserContext    *user = (UserContext*)ctx;
  PetscErrorCode ierr;
  PetscInt       i, j, M, N, xm, ym, xs, ys, num, numi, numj;
  PetscScalar    v[5], Hx, Hy, HydHx, HxdHy;
  MatStencil     row, col[5];
  DM             da;

  PetscFunctionBeginUser;
  ierr  = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  ierr  = DMDAGetInfo(da,0,&M,&N,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

  Hx    = 1.0 / (PetscReal)(M);
  Hy    = 1.0 / (PetscReal)(N);
  HxdHy = Hx/Hy;
  HydHx = Hy/Hx;
  ierr  = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      row.i = i; row.j = j;

      if (i==0 || j==0 || i==M-1 || j==N-1) {
        if (user->bcType == DIRICHLET) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Dirichlet boundary conditions not supported !\n");
        else if (user->bcType == NEUMANN) {
          num=0; numi=0; numj=0;
          if (j!=0) {
            v[num] = -HxdHy;              col[num].i = i;   col[num].j = j-1;
            num++; numj++;
          }
          if (i!=0) {
            v[num] = -HydHx;              col[num].i = i-1; col[num].j = j;
            num++; numi++;
          }
          if (i!=M-1) {
            v[num] = -HydHx;              col[num].i = i+1; col[num].j = j;
            num++; numi++;
          }
          if (j!=N-1) {
            v[num] = -HxdHy;              col[num].i = i;   col[num].j = j+1;
            num++; numj++;
          }
          v[num] = ((PetscReal)(numj)*HxdHy + (PetscReal)(numi)*HydHx); col[num].i = i;   col[num].j = j;
          num++;
          ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);CHKERRQ(ierr);
        }
      } else {
        v[0] = -HxdHy;              col[0].i = i;   col[0].j = j-1;
        v[1] = -HydHx;              col[1].i = i-1; col[1].j = j;
        v[2] = 2.0*(HxdHy + HydHx); col[2].i = i;   col[2].j = j;
        v[3] = -HydHx;              col[3].i = i+1; col[3].j = j;
        v[4] = -HxdHy;              col[4].i = i;   col[4].j = j+1;
        ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (user->bcType == NEUMANN) {
    MatNullSpace nullspace;

    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
    ierr = MatSetNullSpace(J,nullspace);CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
