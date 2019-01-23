/******************************************************************************
 *  FILE NAME  TaskFunctions.cpp
 *  PURPOSE    task functions
 *
 *  SPEC       
 *  NOTES      NONE
 ******************************************************************************/

#include "math.h"

#include "TaskFunctions.h"

/**
 * detailed description
 *
 * @memo    <memo>
 * @return  <return>
 * @param   c
 * @author  Denis
 * @see     nothing
 */
double exact_solution(const coord& c)
{
  return sin(c.rX + c.rY);
}

/**
 * detailed description
 *
 * @memo    <memo>
 * @return  <return>
 * @param   c
 * @author  Denis
 * @see     nothing
 */
double bound_condition(const coord& c)
{
  return exact_solution(c);
}


/**
 * detailed description
 *
 * @memo    <memo>
 * @return  <return>
 * @param   c
 * @author  Denis
 * @see     nothing
 */
double right_part(const coord& c)
{
  return 2.0 * sin(c.rX + c.rY);
}
