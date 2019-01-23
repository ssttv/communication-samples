.PHONY: all

all: project

project: main.cpp Decomposition.cpp FDMApproximation.cpp FDMGrid.cpp Solvers.cpp SparseMatrix.cpp TaskFunctions.cpp
	mpic++ -O3 -Wall -D TEST_CONV -D PRINT_RESULTS main.cpp utils.hpp Decomposition.cpp Decomposition.h FDMApproximation.cpp FDMApproximation.h FDMGrid.cpp FDMGrid.h Solvers.cpp Solvers.h SparseMatrix.cpp SparseMatrix.h TaskFunctions.cpp TaskFunctions.h -o project.out -lm