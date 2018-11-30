/******************************************************************************
 *  FILE NAME  Solvers.cpp
 *  PURPOSE    PCG Method implementation
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#include "stdio.h"
#include "Solvers.h"

int PCGM(const sparse_matrix& mA, vector& vSol, const vector& vF, double rEps, precond_type& prec_data)
{
  double  rAlpha, rBetta, rBRR1, rBRR2;
  size_t  nK, nSize;

  nSize = vF.size();

  vector vR(vF), vP(nSize), vTemp(nSize);

  vSol = 0.0;
  // calc residual
  vTemp.matrix_multiply(mA, vSol);
  vR.add_scale_vector(vTemp, -1.0);

  prec_data.precond_func(vR, vP);
  rBRR1 = vP * vR;
  vTemp = vP;

  nK   = 0;
  while((nK < nSize) && (vTemp.norm_inf() > rEps))
  {
    vTemp.matrix_multiply(mA, vP);
    rAlpha = rBRR1 / (vTemp * vP);
    
    vSol.add_scale_vector(vP, rAlpha);
    vR.add_scale_vector(vTemp, - rAlpha);

    prec_data.precond_func(vR, vTemp);

    rBRR2  = vTemp * vR;
    rBetta = rBRR2 / rBRR1;
    rBRR1  = rBRR2;

    vP *= rBetta;
    vP += vTemp;

    nK ++;
  }

  return int(nK);
}
