
#include "Decomposition.h"

square_decomp::square_decomp(int nDom, int nDomNum, int nN)
{
  int nI = 0, nK = 0, i;
  int nSqrtDomNum = int(sqrt(nDomNum));
  int nSendNum = int(nDom % nSqrtDomNum < nSqrtDomNum - 1) + \
    int(nDom % nSqrtDomNum > 0) + int(nDom / nSqrtDomNum > 0) + \
    int(nDom / nSqrtDomNum < nSqrtDomNum - 1);

  // initialize vectors
  m_vSendInd = std::vector<int>(nSendNum*nN);
  m_vShift   = std::vector<int>(nSendNum);
  m_vRank    = std::vector<int>(nSendNum);
  
  // up
  if (nDom % nSqrtDomNum < nSqrtDomNum - 1)
  {
    m_vRank[nK] = nDom + 1;
    m_vShift[nK] = nI;
    nK++;
    for (i = 0; i < nN; ++i) 
    {
      m_vSendInd[nI] = nN * (i + 1) - 1;
      nI++;
    }
  }
  // down
  if (nDom % nSqrtDomNum > 0)
  {
    m_vRank[nK] = nDom  - 1;
    m_vShift[nK] = nI;
    nK++;
    for (i = 0; i < nN; ++i) 
    {
      m_vSendInd[nI] = nN * i;
      nI++;
    }
    m_vShift[nK] = 0;
  }
  // left
  if (nDom / nSqrtDomNum > 0)
  {
    m_vRank[nK] = nDom - nSqrtDomNum;
    m_vShift[nK] = nI;
    nK++;
    for (i = 0; i < nN; ++i) 
    {
      m_vSendInd[nI] = i;
      nI++;
    }
  }
  // right
  if (nDom / nSqrtDomNum < nSqrtDomNum - 1)
  {
    m_vRank[nK] = nDom + nSqrtDomNum;
    m_vShift[nK] = nI;
    nK++;
    for (i = 0; i < nN; ++i) 
    {
      m_vSendInd[nI] = nN * (nN - 1) + i;
      nI++;
    }
  }
}
