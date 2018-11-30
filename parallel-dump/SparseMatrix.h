/******************************************************************************
 *  FILE NAME  SparseMatrix.h
 *  PURPOSE    sparse matrix
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#pragma once

// ****************************************************************************
// Includes
// ****************************************************************************

#include <valarray>


// ****************************************************************************
// data types

  

/**
 * @memo    allocate matrix
 * @return  a pointer to a matrix
 */
template <typename VALUE_TYPE>
VALUE_TYPE** new_matrix(int size1, int size2)
{
  int i;

  VALUE_TYPE** ppMatrix = new VALUE_TYPE *[size1];
  ppMatrix[0] = new VALUE_TYPE[size1 * size2];

  for (i = 1; i < size1; i ++)
    ppMatrix[i] = ppMatrix[i - 1] + size2;

  return ppMatrix;
}

/**
 * @memo    free a matrix
 */
template <typename VALUE_TYPE>
void delete_matrix(VALUE_TYPE** ppMatrix)
{
  delete[] ppMatrix[0];
  delete[] ppMatrix;
}


class vector;
/**
 * detailed description
 *
 * @memo    <memo>
 * @author  Denis
 * @see     nothing
 */
class sparse_matrix
{
  double **m_mVal;
  int    **m_mInd;
  
  int  m_nSize1;
  int  m_nSize2;
  int  m_nWidth;

  void null_pointers()
  {
    m_mVal = NULL;
    m_mInd = NULL;
  }

  sparse_matrix(const sparse_matrix& matr);
  sparse_matrix& operator=(const sparse_matrix&);

public:
  enum { not_index = -1 };

  sparse_matrix() { null_pointers(); }
  sparse_matrix(int nSize1, int nWidth)
  {
    null_pointers();
    resize(nSize1, nSize1, nWidth);
  }
  sparse_matrix(int nSize1, int nSize2, int nWidth)
  {
    null_pointers();
    resize(nSize1, nSize2, nWidth);
  }
  ~sparse_matrix()
  {
    delete_sm();
  }

  void resize(int nSize1, int nSize2, int nWidth);
  void delete_sm();

  int size1() const { return m_nSize1; }
  int size2() const { return m_nSize2; }
  int width() const { return m_nWidth; }

  double  operator()(int nI, int nJ) const;
  double& add2(int nI, int nJ);

  void del_row(int nRow);
  
  void vector_multiply(const vector& vIn, vector& vOut) const;

  class iterator
  { 
    sparse_matrix& m_refMatrix;
    int            m_nRow;
    int            m_nNonzero;


  public:
    iterator(sparse_matrix& refMatrix) : m_refMatrix(refMatrix), m_nRow(), m_nNonzero() { }
    
    // access interface
    double operator*() const
    {
      return m_refMatrix.m_mVal[m_nRow][m_nNonzero];
    }
    double& operator*()
    {
      return m_refMatrix.m_mVal[m_nRow][m_nNonzero];
    }
    int first_index() const
    {
      return m_nRow;
    }
    int second_index() const
    {
      return m_refMatrix.m_mInd[m_nRow][m_nNonzero];
    }

    // iterator interface
    void First1()
    {
      m_nRow = 0;
      return;
    }
    bool IsDone1() const
    {
      return m_nRow >= m_refMatrix.size1();
    }
    void Next1()
    {
      ++m_nRow;
    }

    void First2()
    {
      m_nNonzero = 0;
    }
    void First2(int nRow)
    {
      m_nRow     = nRow;
      m_nNonzero = 0;
    }
    bool IsDone2() const
    {
      return m_nNonzero >= m_refMatrix.width() ||
             sparse_matrix::not_index == m_refMatrix.m_mInd[m_nRow][m_nNonzero];
    }
    void Next2()
    {
      ++m_nNonzero;
    }
  };
}; // end of class sparse_matrix declaration



/**
 * detailed description
 *
 * @memo    <memo>
 * @author  Denis
 * @see     nothing
 */
class vector : public std::valarray<double>
{
  typedef std::valarray<double>    base_type;

public:
  explicit vector(size_t nSize = 0) : base_type(nSize) { }
  vector(size_t nSize, double rDefault) : base_type(rDefault, nSize) { }

  double  operator*( const vector& vVec ) const;
  vector& operator*=( double rScale );
  vector& operator=( double rRight ) { base_type::operator=(rRight); return *this;}

  double  norm_2()
  {
    return sqrt( (*this) * (*this) );
  }
  double norm_inf()
  {
    double rMax;
    size_t i;

    for (i = 0, rMax = 0.0; i < size(); i++)
    {
      if (rMax < fabs((*this)[i]))
      {
        rMax = fabs((*this)[i]);
      }
    }
    return rMax;
  }

  void    add_scale_vector(const vector& vVec, double rScale);

  void matrix_multiply(const sparse_matrix& mA, vector& vIn)
  {
    mA.vector_multiply(vIn, *this);
  }
}; // end of data_vector




// ****************************************************************************
// include inline functions



