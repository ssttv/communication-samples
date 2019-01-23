#ifndef UTILS
#define UTILS

#include <iostream>

#define MAX_PROC_NUM 1024

template <typename t>
void PrintVect(t* v, size_t s)
{
  for (int i = 0; i < s; ++i)
    std::cout << v[i] << std::endl;
}

static inline double _max(double rX, double rY) 
{
  if (rX > rY)
  {
    return rX;
  }
  else return rY;
}

bool IsSquare(int n)
{
  for (int m = 1; m < MAX_PROC_NUM + 1; ++m)
  {
    if (m*m == n)
    {
      return true;
    }
  }
  return false;
}

#endif