/******************************************************************************
 *  FILE NAME  FDMGrid.h
 *  PURPOSE    grid for finite difference method
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/
#pragma once

#include <vector>

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
  std::vector<coord>  m_vCoords;
  std::vector<bool>   m_vBoundFlag;

  int     m_nN, m_nSize;
  double  m_rH;

public:
  fdm_grid(int nParam);
  ~fdm_grid() { }

  bool node_is_boundary(int nNode) const
  {
    return m_vBoundFlag[nNode];
  }
  void get_neighbours(int nNode, std::vector<int>& vNeighbours) const;
};
