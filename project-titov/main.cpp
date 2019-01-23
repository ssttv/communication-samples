/******************************************************************************
* Solution of the Laplace equation
* Denis Smirnov
* Fedor Garbuzov
*******************************************************************************/


#include <stdio.h>

#include "utils.hpp"
#include "TaskFunctions.h"
#include "FDMApproximation.h"
#include "Solvers.h"

#define MPI_ROOT_PROC 0

int main(int argc, char** argv)
{
  int nGridPart, nBWidth, nNode, nIter, nLocIter;
  double rMaxError;
  int i, j, k, l;
  int nNumProcs, rank;
  int printRank = -1;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nNumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (!IsSquare(nNumProcs))
  {
    if (rank == MPI_ROOT_PROC)
    {
      printf("Error: number of processors must be square number ");
      printf("less than %d.\n", MAX_PROC_NUM);
    }
    MPI_Finalize();
    exit(1);
  }

  // set the global partition number of the square side
  nGridPart = 101; // has to be odd!
  // set maximal number of nonzero matrix elements in row
  nBWidth = 5;

  // create local grid
  fdm_grid grid(nGridPart, nNumProcs, rank);

  // allocate memory for local SLAE matrix, solution and right part vector
  sparse_matrix mA(grid.m_nSize, nBWidth);
  vector vRP(grid.m_nSize), vSol(grid.m_nSize);

  // calculate stiffness matrix and right part
  fdm_slau_assembling(grid, mA, vRP);

  // create decomposition
  square_decomp decomp(rank, nNumProcs, grid.m_nN);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  // solve SLAE
  precond_type prec_data;
  nLocIter = PCGM(mA, vSol, vRP, 1e-9, prec_data, decomp, nNumProcs, rank);
  MPI_Reduce(&nLocIter, &nIter, 1, MPI_INT, MPI_MAX, MPI_ROOT_PROC, MPI_COMM_WORLD);

  // prepare buffers of standart type to send all solutions to the root process
  double *rvGlobSol = new double[grid.m_nGlobN*grid.m_nGlobN];
  double *rvSol     = new double[vSol.size()];
  for (i = 0; i < vSol.size(); ++i)
  {
    rvSol[i] = vSol[i];
  }
  MPI_Gather(rvSol, vSol.size(), MPI_DOUBLE, rvGlobSol, vSol.size(), MPI_DOUBLE,\
      MPI_ROOT_PROC, MPI_COMM_WORLD);

  // postprocessing on the root process
  if (rank == MPI_ROOT_PROC)
  {
    // additional constants
    int nSqrtNumProcs = int(sqrt(nNumProcs));
    int nLocGridPoint = grid.m_nN, nGridPoint = grid.m_nGlobN;
    // solution
    vector vGlobSol(grid.m_nGlobN*grid.m_nGlobN);
    // vector with correct indices
    int *vGlobInd = new int[vGlobSol.size()];
    int nRowInd, nRowIndNew;
    for (i = 0; i < nNumProcs; ++i)
    {
      for (k = 0; k < nLocGridPoint; ++k)
      {
        for (l = 0; l < nLocGridPoint; ++l)
        {
          nRowInd = i/nSqrtNumProcs*nGridPoint*nLocGridPoint + \
              (i%nSqrtNumProcs)*nLocGridPoint + k*nGridPoint + l;
          nRowIndNew = i*nLocGridPoint*nLocGridPoint + k*nLocGridPoint + l;
          vGlobInd[nRowInd] = nRowIndNew;
        }
      }
    }

    // reorder solution
    for (i = 0; i < vGlobSol.size(); ++i)
    {
      vGlobSol[i] = rvGlobSol[vGlobInd[i]];
      //printf("%g %g\n", rvGlobSol[i], vGlobSol[i]);
    }

    // calc max error
    rMaxError = 0.0;
    fdm_grid glob_grid(nGridPart);
    for (nNode = 0; nNode < glob_grid.m_nSize; ++nNode)
    {
      rMaxError = _max(fabs(vGlobSol[nNode] \
        - exact_solution(glob_grid.m_vCoords[nNode])), rMaxError);
    }

    printf("h     = %1.5f\n", grid.m_rH);
    printf("Error = %1.10f\n", rMaxError);
    printf("Iter  = %d\n", nIter);

    delete[] vGlobInd;
  }

  delete[] rvGlobSol;
  delete[] rvSol;

  MPI_Finalize();
  return 0;
}

