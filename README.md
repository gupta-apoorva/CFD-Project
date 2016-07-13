 
                                                            Project Summary (Group 3)

Question: To integrate the Energy transport equation into the NAvier Stokes and solving the system using PETSc instead of SOR.  

						Validation Case 1 (Page 138 , Numerical Simulation in Fluid Dynamics):

In this case the left wall is heated and the right wall is cold. The normalised values are given in the book and so we have implemented them in the code.

How to RUN:  

1> Make the code using "make all"

2> Run the code using one of the options

./sim -da_grid_x 50 -da_grid_y 50 -pc_type mg -pc_mg_levels 1 -mg_levels_0_pc_type lu -mg_levels_0_pc_factor_shift_type NONZERO -ksp_monitor 

This solves the poissons part of the problem using Multigrid with the last level being solved using LU decomposition.  -da_grid_x and -da_grid_y are the grid sizes of the problem in the respected direction.


						
						Validation Case 2 (Page 143 , Numerical Simulation in Fluid Dynamics):


In this case the bottom wall is heated and the top wall is cold. The normalised values are given in the book and we have implemened those.

How to RUN:

1> Make the code using "make all"

2> Run the code using one of the options

    ./sim -da_grid_x 65 -da_grid_y 17 -pc_type mg -pc_mg_levels 1 -mg_levels_0_pc_type lu -mg_levels_0_pc_factor_shift_type NONZERO -ksp_monitor

This solves the poissons part of the problem using Multigrid with the last level being solved using LU decomposition.  -da_grid_x and -da_grid_y are the grid sizes of the problem in the respected direction

    ./sim -da_grid_x 65 -da_grid_y 17 -pc_type mg -da_refine 0 -ksp_atol 0.0001 -ksp_monitor

This solves the poissons part of the problem using Multigrid with the refinement level being ZERO .  -da_grid_x and -da_grid_y are the grid sizes of the problem in the respected direction on the coarsest level and -ksp_atol is the absolute tolerance we want to give to the solver.

    ./sim -da_grid_x 33 -da_grid_y 9 -pc_type mg -da_refine 1 -ksp_atol 0.0001 -ksp_monitor

This solves the poissons part of the problem using Multigrid with the refinement level being ONE .  -da_grid_x and -da_grid_y are the grid sizes of the problem in the respected direction on the coarsest level and -ksp_atol is the absolute tolerance we want to give to the solver.

    ./sim -da_grid_x 17 -da_grid_y 5 -pc_type mg -da_refine 2 -ksp_atol 0.0001 -ksp_monitor

This solves the poissons part of the problem using Multigrid with the refinement level being TWO .  -da_grid_x and -da_grid_y are the grid sizes of the problem in the respected direction on the coarsest level and -ksp_atol is the absolute tolerance we want to give to the solver.

This solves the poissons part of the problem using Multigrid with the refinement level being THREE .  -da_grid_x and -da_grid_y are the grid sizes of the problem in the respected direction on the coarsest level and -ksp_atol is the absolute tolerance we want to give to the solver.

    ./sim -da_grid_x 9 -da_grid_y 3 -pc_type mg -da_refine 3 -ksp_atol 0.0001 -ksp_monitor


