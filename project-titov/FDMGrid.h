/******************************************************************************
 *  FILE NAME  FDMGrid.h
 *  PURPOSE    grid for finite difference method
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/
#ifndef FDMGRID
#define FDMGRID

#include <vector>
#include "math.h"

struct coord
{
   double  rX, rY;
};

/**
 * detailed description
 *
 * @memo    <memo>
 * @author  Denis
 * @see     nothing
 */
class fdm_grid
{
public:
  std::vector<coord> m_vCoords;
  std::vector<bool>  m_vBoundFlag;
  std::vector<bool>  m_vInterDomUpFlag;
  std::vector<bool>  m_vInterDomDownFlag;
  std::vector<bool>  m_vInterDomLeftFlag;
  std::vector<bool>  m_vInterDomRightFlag;

  int     m_nN, m_nSize, m_nGlobN, m_nDom, m_nDomNum;
  double  m_rH;

public:
  //fdm_grid(int nParam);
  fdm_grid(int nParam, int nDomNum = 1, int nDom = 0);
  ~fdm_grid() { }

  bool node_is_boundary(int nNode) const
  {
    if (nNode < m_nSize)  return m_vBoundFlag[nNode];
    return false;
  }
  void get_neighbours(int nNode, std::vector<int>& vNeighbours) const;
};

#endif