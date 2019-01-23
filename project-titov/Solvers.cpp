/******************************************************************************
 *  FILE NAME  Solvers.cpp
 *  PURPOSE    PCG Method implementation
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#include "stdio.h"
#include "Solvers.h"


int PCGM(const sparse_matrix& mA, vector& vSol, const vector& vF, double rEps, \
  precond_type& prec_data, decomposition& decomp, int nNumProcs, int rank)
{
  double rAlpha, rBeta, rBRR1, rLocBRR1, rBRR2, rLocBRR2, rNorm, rLocNorm;
  double rTmp, rLocTmp;
  size_t nK, nSize = vF.size(), nExtSize = decomp.m_vRank.size();
  int i, j, k, len, nExtVectLen = nSize;

  MPI_Status  vStatus[nExtSize];
  MPI_Request vSendReqArr[nExtSize];
  MPI_Request vRecvReqArr[nExtSize];

  // buffers to send elements and to store received elements
  double **vvSendVal = (double**)malloc(sizeof(double*)*nExtSize);
  double **vvRecvVal = (double**)malloc(sizeof(double*)*nExtSize);
  for (i = 0; i < nExtSize; ++i)
  {
    if (i == nExtSize - 1)
        len = decomp.m_vSendInd.size() - decomp.m_vShift[nExtSize - 1];
      else
        len = decomp.m_vShift[i + 1] - decomp.m_vShift[i];
    vvSendVal[i] = (double*)malloc(sizeof(double)*len);
    vvRecvVal[i] = (double*)malloc(sizeof(double)*len);
    nExtVectLen += len;
  }
 
  // allocate memory for extended vector and other vectors 
  vector vExtVal(nExtVectLen);
  vector vR(vF), vP(nSize), vTemp(nSize);
  vSol  = 0.0;
  vTemp = 0.0;

  prec_data.precond_func(vR, vP);
  rLocBRR1 = vP * vR;
  MPI_Allreduce(&rLocBRR1, &rBRR1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  vTemp = vP;

  // main loop
  nK   = 0;
  rNorm = 10000.0;
  while((nK < nSize*nNumProcs) && (rNorm > rEps))
  {
    // exchange interdomain node values
    for (i = 0; i < nExtSize; ++i)
    {
      if (i == nExtSize - 1)
        len = decomp.m_vSendInd.size() - decomp.m_vShift[nExtSize - 1];
      else
        len = decomp.m_vShift[i + 1] - decomp.m_vShift[i];
      for (j = 0; j < len; ++j)
        vvSendVal[i][j] = vP[decomp.m_vSendInd[decomp.m_vShift[i] + j]];
      MPI_Isend(vvSendVal[i], len, MPI_DOUBLE, decomp.m_vRank[i], \
        decomp.m_vRank[i], MPI_COMM_WORLD, vSendReqArr + i);
      MPI_Irecv(vvRecvVal[i], len, MPI_DOUBLE, decomp.m_vRank[i], \
        rank, MPI_COMM_WORLD, vRecvReqArr + i);
    }

    // wait for all requests
    MPI_Waitall(nExtSize, vSendReqArr, vStatus);
    MPI_Waitall(nExtSize, vRecvReqArr, vStatus);
    MPI_Barrier(MPI_COMM_WORLD);

    // copy vP to the extended vector
    for (k = 0; k < vP.size(); ++k)
      vExtVal[k] = vP[k];

    // add received values to the extended vector
    for (i = 0; i < nExtSize; ++i)
    {
      len = decomp.m_vShift[i + 1] - decomp.m_vShift[i];
      if (i == nExtSize - 1)
        len = decomp.m_vSendInd.size() - decomp.m_vShift[nExtSize - 1];
      for (j = 0; j < len; ++j)
      {
        vExtVal[k] = vvRecvVal[i][j];
        ++k;
      }
    }

    // multiply matrix by the extended vector
    vTemp.matrix_multiply(mA, vExtVal);

    rLocTmp = vTemp * vP;
    MPI_Allreduce(&rLocTmp, &rTmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    rAlpha = rBRR1 / rTmp;
    
    vSol.add_scale_vector(vP, rAlpha);
    vR.add_scale_vector(vTemp, - rAlpha);

    prec_data.precond_func(vR, vTemp);

    rLocBRR2  = vTemp * vR;
    MPI_Allreduce(&rLocBRR2, &rBRR2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    rBeta = rBRR2 / rBRR1;
    rBRR1  = rBRR2;

    vP *= rBeta;
    vP += vTemp;

    nK ++;

    rLocNorm = vTemp.norm_inf();
    MPI_Allreduce(&rLocNorm, &rNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }

  for (i = 0; i < nExtSize; ++i)
  {
    free(vvSendVal[i]);
    free(vvRecvVal[i]);
  }
  free(vvSendVal);
  free(vvRecvVal);

  return int(nK);
}
