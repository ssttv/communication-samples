/******************************************************************************
 *  FILE NAME  Solvers.h
 *  PURPOSE    PCG Method
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#pragma once

#include "Preconditioners.h"

extern int PCGM(const sparse_matrix& mA, vector& vSol, const vector&vF, double rEps, precond_type& prec_data);