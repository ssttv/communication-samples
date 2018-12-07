/******************************************************************************
 *  FILE NAME  SparseMatrix.cpp
 *  PURPOSE    sparse matrix implementation
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

// ****************************************************************************
// Includes
// ****************************************************************************
#include <iostream>

#include "SparseMatrix.h"


// ****************************************************************************
// data types



/**
 * detailed description
 *
 * @memo    <memo>
 * @param   nSize1
 * @param   nSize2
 * @param   nWidth
 * @author  Denis
 * @see     nothing
 */
void sparse_matrix::resize(int nSize1, int nSize2, int nWidth)
{
  int nI, nK;

  delete_sm();

  m_mVal = new_matrix<double>(nSize1, nWidth);
  m_mInd = new_matrix<int>(nSize1, nWidth);

  m_nSize1  = nSize1;
  m_nSize2  = nSize2;
  m_nWidth  = nWidth;

  for (nI = 0; nI < nSize1; nI++)
  {
     for (nK = 0; nK < nWidth; nK++)
     {
       m_mVal[nI][nK] = 0.0;
       m_mInd[nI][nK] = not_index;
     }
  }
}

/**
 * detailed description
 *
 * @memo    <memo>
 * @author  Denis
 * @see     nothing
 */
void sparse_matrix::delete_sm()
{
  if (NULL != m_mVal && NULL != m_mInd)
  {
    delete_matrix(m_mVal);
    delete_matrix(m_mInd);
  }
}

/**
 * detailed description
 *
 * @memo    <memo>
 * @return  <return>
 * @param   
 * @author  Denis
 * @see     nothing
 */
double sparse_matrix::operator()(int nI, int nJ) const
{
  int nK, nSecondIdx;

  for (nK = 0; nK < m_nWidth; nK++)
  {
    nSecondIdx = m_mInd[nI][nK];
    if ( nSecondIdx == nJ )
    {
      return m_mVal[nI][nK];
    }
  }

  return 0;
}

/**
 * detailed description
 *
 * @memo    <memo>
 * @return  <return>
 * @param   
 * @author  Denis
 * @see     nothing
 */
double& sparse_matrix::add2(int nI, int nJ)
{
  int nK, nSecondIdx;

  for (nK = 0; nK < m_nWidth; nK++)
  {
    nSecondIdx = m_mInd[nI][nK];
    if ( nSecondIdx == not_index || nSecondIdx == nJ )
    {
      m_mInd[nI][nK] = nJ;

      return m_mVal[nI][nK];
    }
  }

  // TODO: throw exception
  std::cout << "incorrect width of the matrix" <<std::endl;
  return m_mVal[0][0];
}

/**
 * detailed description
 *
 * @memo    <memo>
 * @return  <return>
 * @param   nRow
 * @author  Denis
 * @see     nothing
 */

void sparse_matrix::del_row(int nRow)
{
  int nK;

  for (nK = 0; nK < m_nWidth; nK++)
  {
    m_mInd[nRow][nK] = not_index;
    m_mVal[nRow][nK] = 0;
  }
}

/**
 * detailed description
 *
 * @memo    <memo>
 * @return  <return>
 * @param   vIn
 * @param   vOut
 * @author  Denis
 * @see     nothing
 */

void sparse_matrix::vector_multiply(const vector& vIn, vector& vOut) const
{
  int nI, nJ, nIdx;

  vOut = 0.0;

  for (nI = 0; nI < m_nSize1; nI++)
  {
    for (nJ = 0; nJ < m_nWidth; nJ++)
    {
      nIdx = m_mInd[nI][nJ];
      if (not_index == nIdx)
      {
        break;
      }
      vOut[nI] += m_mVal[nI][nJ] * vIn[nIdx];
    }
  }

  return;
}

/**
 * detailed description
 *
 * @memo    dot product
 * @return  none
 * @param   vVec
 * @author  Denis
 * @see     nothing
 */

double vector::operator*( const vector& vVec ) const
{
  double  rDotProduct;
  size_t  nI;

  rDotProduct = 0.0;
  for (nI = 0; nI < size(); ++nI)
  {
    rDotProduct += (*this)[nI] * vVec[nI];
  }

  return rDotProduct;
} // end of operator*

/**
 * detailed description
 *
 * @memo    dot product
 * @return  none
 * @param   vVec
 * @author  Denis
 * @see     nothing
 */
vector& vector::operator*=( double rScale )
{
  size_t  nI;

  for (nI = 0; nI < this->size(); ++nI)
  {
    (*this)[nI] *= rScale;
  }

  return *this;
} // end of operator*=




/**
 * detailed description
 *
 * @memo    this += vVec * rScale
 * @return  none
 * @param   vVec    - vector
 * @param   rScale  - scale
 * @author  Denis
 * @see     nothing
 */

void vector::add_scale_vector( const vector& vVec, double rScale )
{
  size_t  nI;

  for (nI = 0; nI < this->size(); ++nI)
  {
    (*this)[nI] += vVec[nI] * rScale;
  }
} // end of add_scale_vector






