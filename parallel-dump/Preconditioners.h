/******************************************************************************
 *  FILE NAME  Preconditioners.h
 *  PURPOSE    precondition operator
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#pragma once

#include "SparseMatrix.h"

class precond_type
{
  public:
    void precond_func(const vector& vIn, vector& vOut)
    {
      vOut = vIn;
    }
};