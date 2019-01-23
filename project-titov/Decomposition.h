#ifndef DECOMPOSITION
#define DECOMPOSITION

#include <vector>
#include <math.h>

class decomposition
{
public:
  std::vector<int> m_vSendInd;
  std::vector<int> m_vShift;
  std::vector<int> m_vRank;
public:
  decomposition() { }
  ~decomposition() { }
};

class square_decomp : public decomposition
{
public:
  square_decomp(int nDom, int nDomNum, int nN);
  ~square_decomp() {}
};

#endif