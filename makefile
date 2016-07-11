#CFLAGS	         =
#FFLAGS	         =
#CPPFLAGS         =
#FPPFLAGS         =
#LOCDIR           = #src/ksp/ksp/examples/tutorials/
#EXAMPLESC        = multigrid.c
#EXAMPLESF        = 
#MANSEC           = KSP
#CLEANFILES       = rhs.vtk solution.vtk
#NP               = 1
OBJ = 	helper.o\
      	init.o\
	boundary_val.o\
      	uvp.o\
	sor.o\
      	main.o\
      	visual.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

#multigrid: multigrid.o chkopts
	#gcc -o multigrid multigrid.o ${PETSC_KSP_LIB}
	
all: multigrid.o chkopts $(OBJ)
	gcc -o sim multigrid.o $(OBJ)  -lm ${PETSC_KSP_LIB}
	${RM} multigrid.o $(OBJ)


include ${PETSC_DIR}/lib/petsc/conf/test

#-${CLINKER}
