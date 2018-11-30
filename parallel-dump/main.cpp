/******************************************************************************
 *  FILE NAME  FDMGrid.h
 *  PURPOSE    solving of the Laplace's equation
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/


#include <stdio.h>

#include "TaskFunctions.h"
#include "FDMApproximation.h"
#include "Solvers.h"

int main()
{
  int            nGridParam, nBWidth, nNode, nIter;
  double         rMaxError;

  // set the partition number of the square side
  nGridParam = 10;
  // set maximal number of nonzero matrix elements in row
  nBWidth    = 5;

  // create grid
  fdm_grid   grid(nGridParam);

  // allacate SLAE
  // - sparse matrix
  sparse_matrix  mA(grid.m_nSize, nBWidth);
  // - solution and right part
  vector vSol(grid.m_nSize), vRP(grid.m_nSize);

  // calculate stiffness matrix and right part
  fdm_slau_assembling(grid, mA, vRP);

  // solve SLAE
  precond_type prec_data;
  nIter = PCGM(mA, vSol, vRP, 1e-9, prec_data);

  // calc max error
  rMaxError = 0.0;
  for (nNode = 0; nNode < grid.m_nSize; ++nNode)
  {
    rMaxError = std::max(fabs(vSol[nNode] - exact_solution(grid.m_vCoords[nNode])), rMaxError);
  }

  printf("h     = %1.5f\n", grid.m_rH);
  printf("Error = %1.10f\n", rMaxError);
  printf("Iter  = %d\n", nIter);
  
  /* 
  // print nonzero matrix elements using iterator (variant 1)
  sparse_matrix::iterator itA(mA);
  for (itA.First1(); !itA.IsDone1(); itA.Next1())
  {
    printf("\n");
    for (itA.First2(); !itA.IsDone2(); itA.Next2())
    {
      printf("%5.1f(%d %d) ", *itA, itA.first_index(), itA.second_index());
    }
  }*/

  /*
  // print nonzero matrix elements using iterator (variant 2)
  sparse_matrix::iterator itA(mA);
  int    nI;
  for (nI = 0; nI < mA.size1(); nI++)
  {
    printf("\n");
    for (itA.First2(nI); !itA.IsDone2(); itA.Next2())
    {
      printf("%5.1f(%d %d) ", *itA, itA.first_index(), itA.second_index());
    }
  }*/

  return 0;
}

