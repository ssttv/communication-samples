/******************************************************************************
 *  FILE NAME  FDMGrid.cpp
 *  PURPOSE    grid for finite difference method
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#include <iostream>

#include "FDMGrid.h"


fdm_grid::fdm_grid(int nParam, int nDomNum, int nDom)
{
  // NOTE: number of domains (nDomNum) has to be a square number
  int     nXLine, nI, nJ, nDomCol, nDomRow;
  int     nSqrtDomNum = int(sqrt(nDomNum));
  double  rXCoord, rMinXCoord, rMinYCoord;
  
  nDomCol = nDom / nSqrtDomNum;
  nDomRow = nDom % nSqrtDomNum;

  m_rH    = 1.0 / nParam;
  m_nGlobN = nParam + 1;
  m_nN    = m_nGlobN/nSqrtDomNum;
  m_nSize = m_nN * m_nN;
  m_nDom = nDom;
  m_nDomNum = nDomNum;

  rMinXCoord = nDomCol*m_rH*m_nN;
  rMinYCoord = nDomRow*m_rH*m_nN;

  m_vCoords.resize(m_nSize);
  m_vBoundFlag.resize(m_nSize, false);
  m_vInterDomUpFlag.resize(m_nSize, false);
  m_vInterDomDownFlag.resize(m_nSize, false);
  m_vInterDomLeftFlag.resize(m_nSize, false);
  m_vInterDomRightFlag.resize(m_nSize, false);

  for (nI = 0; nI < m_nN; ++nI)
  {
    // set coordinates
    nXLine  = nI * m_nN;
    rXCoord = nI * m_rH + rMinXCoord;
    for (nJ = 0; nJ < m_nN; ++nJ)
    {
      m_vCoords[nXLine + nJ].rX = rXCoord;
      m_vCoords[nXLine + nJ].rY = nJ * m_rH + rMinYCoord;
    }

    // mark boundary nodes
    // - horizontal boundary
    if (nDomRow == 0)
      m_vBoundFlag[ nXLine             ] = true;
    if (nDomRow == nSqrtDomNum - 1)
      m_vBoundFlag[ nXLine  + m_nN - 1 ] = true;
    // - vertical boundary
    if (nDomCol == 0)
      m_vBoundFlag[ 0                  + nI ] = true;
    if (nDomCol == nSqrtDomNum - 1)
      m_vBoundFlag[ (m_nN - 1) * m_nN  + nI ] = true;

    // mark interdomain boundary nodes
    m_vInterDomUpFlag[nXLine + m_nN - 1] = (nDomRow < nSqrtDomNum - 1);
    m_vInterDomDownFlag[nXLine] = (nDomRow > 0);
    m_vInterDomLeftFlag[0 + nI] = (nDomCol > 0);
    m_vInterDomRightFlag[(m_nN - 1) * m_nN  + nI] = (nDomCol < nSqrtDomNum - 1);
  }
}



void fdm_grid::get_neighbours(int nNode, std::vector<int>& vNeighbours) const
{
  int nExtNum = 0;
  int nSqrtDomNum = int(sqrt(m_nDomNum));
  int nDomCol = m_nDom / nSqrtDomNum;
  int nDomRow = m_nDom % nSqrtDomNum;

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
    return;
  }

  if (m_vInterDomUpFlag[nNode])
  {
    vNeighbours[1] = m_nSize + nNode/m_nN;
  }
  if (m_vInterDomDownFlag[nNode])
  {
    nExtNum = int(nDomRow < nSqrtDomNum - 1);
    vNeighbours[0] = m_nSize + nExtNum*m_nN + nNode/m_nN;
  }
  if (m_vInterDomLeftFlag[nNode])
  {
    nExtNum = int(nDomRow < nSqrtDomNum - 1) + int(nDomRow > 0);
    vNeighbours[2] = m_nSize + nExtNum*m_nN + nNode%m_nN;
  }
  if (m_vInterDomRightFlag[nNode])
  {
    nExtNum = int(nDomRow < nSqrtDomNum - 1) + int(nDomRow > 0) + int(nDomCol > 0);
    vNeighbours[3] = m_nSize + nExtNum*m_nN + nNode%m_nN;
  }
}

