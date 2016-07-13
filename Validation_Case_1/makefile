OBJ = 	helper.o\
      	init.o\
		boundary_val.o\
      	uvp.o\
		sor.o\
      	main.o\
      	visual.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
	
all: multigrid.o chkopts $(OBJ)
	gcc -o sim multigrid.o $(OBJ)  -lm ${PETSC_KSP_LIB} -fopenmp -Wall -pg 


include ${PETSC_DIR}/lib/petsc/conf/test

#-${CLINKER}
