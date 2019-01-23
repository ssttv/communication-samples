/******************************************************************************
 *  FILE NAME  TaskFunctions.h
 *  PURPOSE    task functions
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/
#ifndef TASKFUNC
#define TASKFUNC

#include "FDMGrid.h"

extern double exact_solution(const coord& c);
extern double bound_condition(const coord& c);
extern double right_part(const coord& c);

#endif
