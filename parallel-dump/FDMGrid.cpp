/******************************************************************************
 *  FILE NAME  FDMGrid.cpp
 *  PURPOSE    grid for finite difference method
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#include <iostream>

#include "FDMGrid.h"

fdm_grid::fdm_grid(int nParam)
{
  int     nXLine, nI, nJ;
  double  rXCoord;

  m_rH    = 1.0 / nParam;
  m_nN    = nParam + 1;
  m_nSize = m_nN * m_nN;

  m_vCoords.resize(m_nSize);
  m_vBoundFlag.resize(m_nSize, false);


  for (nI = 0; nI < m_nN; ++nI)
  {
    // set coordinates
    nXLine  = nI * m_nN;
    rXCoord = nI * m_rH;
    for (nJ = 0; nJ < m_nN; ++nJ)
    {
      m_vCoords[nXLine + nJ].rX = rXCoord;
      m_vCoords[nXLine + nJ].rY = nJ * m_rH;
    }

    // mark boundary nodes
    // - horizontal boundary
    m_vBoundFlag[ nXLine             ] = true;
    m_vBoundFlag[ nXLine  + m_nN - 1 ] = true;
    // - vertical boundary
    m_vBoundFlag[ 0                  + nI ] = true;
    m_vBoundFlag[ (m_nN - 1) * m_nN  + nI ] = true;
  }
}



void fdm_grid::get_neighbours(int nNode, std::vector<int>& vNeighbours) const
{
 if (4 == vNeighbours.size())
 {
   vNeighbours[0] = nNode - 1;
   vNeighbours[1] = nNode + 1;
   vNeighbours[2] = nNode - m_nN;
   vNeighbours[3] = nNode + m_nN;
 }
 else
 {
   std::cout << "invalide neieghbours vector size" << std::endl;
 }
}

