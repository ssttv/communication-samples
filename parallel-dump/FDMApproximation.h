/******************************************************************************
 *  FILE NAME  FDMApproximation.h
 *  PURPOSE    Finite difference method
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/
#pragma once

#include "FDMGrid.h"
#include "SparseMatrix.h"

extern void fdm_slau_assembling(const fdm_grid& grid, sparse_matrix& mSM, vector& vRP);
