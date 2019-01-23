/******************************************************************************
 *  FILE NAME  FDMApproximation.cpp
 *  PURPOSE    Finite difference method
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#include "TaskFunctions.h"

#include "FDMApproximation.h"


void fdm_slau_assembling(const fdm_grid& grid, sparse_matrix& mSM, vector& vRP)
{
  double  rH2;
  int     nNode, nI;
  
  std::vector<int> vNeighbours(4);

  rH2 = grid.m_rH * grid.m_rH;

  // assemble SLAE 
  for (nNode = 0; nNode < grid.m_nSize; ++nNode)
  {
    if (grid.node_is_boundary(nNode))
    {
      // write boundary condition
      vRP[nNode] = bound_condition(grid.m_vCoords[nNode]);
      mSM.add2(nNode, nNode) = 1.0;
    }
    else
    {
      // write finite difference equation
      vRP[nNode] = rH2 * right_part(grid.m_vCoords[nNode]);
      mSM.add2(nNode, nNode) = 4.0;
      
      grid.get_neighbours(nNode, vNeighbours);
      for (nI = 0; nI < 4; ++nI)
      {
        if (grid.node_is_boundary(vNeighbours[nI]))
        {
          // take into account boundary condition for SLAE symmetrization
          vRP[nNode] -= -1.0 * bound_condition(grid.m_vCoords[vNeighbours[nI]]);
        }
        else
        {
          mSM.add2(nNode, vNeighbours[nI]) = -1.0;
        }
      }
    }
  }
} // end of fdm_slau_assembling

