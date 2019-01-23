/******************************************************************************
 *  FILE NAME  Solvers.h
 *  PURPOSE    PCG Method
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#ifndef SOLVERS
#define SOLVERS

#include "Preconditioners.h"
#include "Decomposition.h"

extern int PCGM(const sparse_matrix& mA, vector& vSol, const vector&vF, \
    double rEps, precond_type& prec_data, decomposition& decomp, int nNumProcs,\
    int rank);

#endif