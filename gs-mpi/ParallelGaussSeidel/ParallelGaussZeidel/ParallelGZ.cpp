#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>


// computing accuracy
#define EPS 0.1
int ProcNum = 0;      // Number of available processes 
int ProcRank = -1;     // Rank of current process

// Function for formatted matrix output
void PrintMatrix (double* pMatrix, int RowCount, int ColCount) {
  int i, j; // Loop variables
  for (i=0; i<RowCount; i++) {
    for (j=0; j<ColCount; j++)
      printf("%7.4f ", pMatrix[i*ColCount+j]);
    printf("\n");
  }
}

void TestDistribution(double* pMatrix, double* pProcRows,
    int Size, int RowNum) {
  if (ProcRank == 0) {
    printf("Initial Matrix: \n");
    PrintMatrix(pMatrix, Size, Size);    
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i=0; i<ProcNum; i++) {
    if (ProcRank == i) {
      printf("\nProcRank = %d \n", ProcRank);
     // fprintf(" Matrix Stripe:\n");
      PrintMatrix(pProcRows, RowNum, Size);      
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
// Function for distribution of the initial objects between the processes
void DataDistribution(double* pMatrix, double* pProcRows, int Size,int RowNum) {
 MPI_Scatter(pMatrix, RowNum*Size, MPI_DOUBLE, pProcRows, RowNum*Size,MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// Function for Gauss_Zeidel algoritm
void ResultCalculation(double* pMatrix, int Size, int &Step) {
  int i, j;  // Loop variables  
  double temp;
  double dmax;
  double dm; 
  Step = 0;   
  do
  {
    dmax = 0;
    for (i = 1; i < Size - 1; i++)
      for(j = 1; j < Size - 1; j++)
      {
        temp = pMatrix[Size * i + j];
        pMatrix[Size * i + j] = 0.25 * (pMatrix[Size * i + j + 1] + pMatrix[Size * i + j - 1] + 
          pMatrix[Size * (i + 1) + j] + pMatrix[Size * (i - 1) + j]);
        dm = fabs(pMatrix[Size * i + j] - temp);
        if (dmax < dm)
          dmax = dm;        
      }    
      Step++;
  }
  while (dmax > EPS);  
}

// Function for computational process termination
void ProcessTermination(double* pMatrix, double* pProcRows) {
  if (ProcRank == 0)
	delete [] pMatrix;
  delete [] pProcRows;
}

void DummyDataInitialization (double* pMatrix, int Size) {
  int i, j;  // Loop variables
  double h = 1.0 / (Size - 1);
  for (i=0; i<Size; i++) {
    for (j=0; j<Size; j++)
      pMatrix[i*Size+j] = 0;
  }
  for (i = 0; i < Size; i++)
  {
    pMatrix[Size * i + 0] = 100 - 200 * i * h;
    pMatrix[Size * i + Size-1] = -100 + 200 * i * h;
    pMatrix[Size * 0 + i] = 100 - 200 * i * h;
    pMatrix[Size * (Size-1) + i] = -100 + 200 * i * h;
  }
}
// Function for memory allocation and definition of object’s elements
void ProcessInitialization (double* &pMatrix, double* &pProcRows, int & Size, int &RowNum) {    
  if (ProcRank == 0) {
	  fprintf(stdout, "\nEnter the size of initial objects: ");
	  fflush(stdout);
	  scanf_s("%d", &Size);
  }
  MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (Size < ProcNum) {
    if (ProcRank == 0) {
      fprintf(stdout, "Chosen size of the objects = %d \n", Size); 
      fprintf(stdout, "Size of the objects must be greater than number of processes! \n");
      fflush(stdout);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
    return;
  }
  if (Size%ProcNum != 0) {
    if (ProcRank == 0) {
      fprintf(stdout, "Chosen size of the objects = %d \n", Size);         
      fprintf(stdout, "Number of processes must be multiple of size of objects! \n");
      fflush(stdout);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
    return;
  }
  // Define the number of matrix rows stored on each process
  RowNum = (Size-2)/(ProcNum-1)+2;
  pProcRows = new double [RowNum*Size];
  // Define the values of initial objects’ elements
  if (ProcRank == 0) {
    // Initial matrix exists only on the pivot process
    pMatrix = new double [Size*Size];
    // Values of elements are defined only on the pivot process
    DummyDataInitialization(pMatrix, Size);
  }

}
// Function for simple definition of matrix and vector elements

void main(int argc, char* argv[]) {
  double* pMatrix;  // The first argument - initial matrix
  int Size;		    // Sizes of initial matrix and vector
  int Step;
  int RowNum;          // Number of rows in matrix stripe
  double* pProcRows;   // Stripe of the matrix on current process
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
  MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

  if(ProcRank == 0) {
    printf("Parallel Gauss - Zeidel program\n");
    fflush(stdout);
	// Memory allocation and data initialization
    ProcessInitialization(pMatrix, pProcRows, Size, RowNum);
  }
  // Distributing the initial objects between the processes
  DataDistribution(pMatrix, pProcRows, Size, RowNum);
  // Distribution test
  TestDistribution(pMatrix, pProcRows, Size, RowNum);

  if (ProcRank == 0)
    printf ("Number of abailable processes = %d\n", ProcNum);
    printf ("Rank of current process = %d\n", ProcRank);
		
  ProcessTermination(pMatrix, pProcRows);
  MPI_Finalize();

}