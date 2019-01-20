#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <mpi.h>
//#include <libconfig.h>

/**
 * @memo    allocate matrix
 * @return  a pointer to a matrix
 */
template <typename VALUE_TYPE>
VALUE_TYPE** malloc_matrix(int size1, int size2)
{
  int i;

  VALUE_TYPE** ppMatrix = (VALUE_TYPE**)malloc(size1 * sizeof(VALUE_TYPE*));
  ppMatrix[0] = (VALUE_TYPE*)malloc(size1 * size2 * sizeof(VALUE_TYPE));

  for (i = 1; i < size1; i ++)
    ppMatrix[i] = ppMatrix[i - 1] + size2;

  return ppMatrix;
}

/**
 * @memo    free a matrix
 */
template <typename VALUE_TYPE>
void free_matrix(VALUE_TYPE** ppMatrix)
{
  free(ppMatrix[0]);
  free(ppMatrix);
}

/**
 * @memo    return pointer to data array of matrix
 * @return  pointer to matrix
 */
template <typename VALUE_TYPE>
VALUE_TYPE* point_to_data(VALUE_TYPE** ppMatrix)
{
  if (ppMatrix)
  {
    return ppMatrix[0];
  }
  else
  {
    return NULL;
  }
}

int main(int agrc, char* argv[])
{
  double *x;
  int ProcNum; 
  int ProcRank; 
  int i, j, sum, RestSize;
  int n;
	double eps;
  double **A = NULL;
  double **ALocal = NULL;
  double *f = NULL;
  double *fLocal = NULL;
  double residual, LocalResidualNorm, ResidualNorm;
  int *sendCount;
  int *displs;

  // для вывода на экран без буферизации
  setvbuf(stdout, 0, _IONBF, 0);


	MPI_Init(&agrc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
  {
		printf("Input SLAE dimension = ");
		if (scanf("%d", &n) > 0)
      printf("OK\n");
		printf("Input eps = ");
		if (scanf("%lf", &eps) > 0)
      printf("OK\n");
		printf("\n");
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  // Размерности даннных здесь и далее я изменил с n+1 на n
	x = (double*)malloc(n*sizeof(double));
	for (i = 0; i < n; i++)
  {
		x[i] = 0;
  }

  // Глобальную матрицу собираем только на нулевом процессе.
  // Я не понял, по какой логике собиралась матрица в вашем коде, и поэтому сделал свою простую сборку СЛАУ,
  // решение которой - вектор из едениц.
  if (ProcRank == 0)
  {
    // выделяем память под глобальную матрицу и правую часть
    A = malloc_matrix<double>(n, n);
    f = (double*)malloc(n*sizeof(double));

    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        A[i][j] = -1;
      }
      A[i][i] = n;
      f[i] = 1;
    }
  } 

  sendCount = (int*)malloc(sizeof(int)*ProcNum);
  displs = (int*)malloc(sizeof(int)*ProcNum);
  
  // Задаем вспомогательные вектора для функций MPI_Scatterv и MPI_Allgatherv.
  // Добавил RestSize для случая, когда n не делится нацело на ProcNum.
  sum = 0;
  RestSize = n;
  for (int i = 0; i < ProcNum; i++) {
      sendCount[i] = RestSize / (ProcNum - i);
      displs[i]    = sum;

      RestSize -= sendCount[i];
      sum      += sendCount[i];
  }
  // создаем тип сьтрока для раздачи матрицы процессам
  MPI_Datatype   ROW;
  MPI_Type_contiguous(n, MPI_DOUBLE, &ROW);
  MPI_Type_commit(&ROW);

  // выделяем память под локальную матрицы и правую часть
  ALocal = malloc_matrix<double>(sendCount[ProcRank], n);
  fLocal = (double*)malloc(sendCount[ProcRank]*sizeof(double));

  // Раздаем
  // - матрицу раздаем строчками с помощью типа ROW
  MPI_Scatterv(point_to_data(A), sendCount, displs, ROW, point_to_data(ALocal), sendCount[ProcRank], ROW, 0, MPI_COMM_WORLD);
  // - правую часть раздаем с типом MPI_DOUBLE
  MPI_Scatterv(f, sendCount, displs, MPI_DOUBLE, fLocal, sendCount[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	   
  if (ProcRank == 0) printf("Computation in progress...\n\n");

  // Метод Зейделя
  while(1)
  {
    LocalResidualNorm = 0;
    for (i = 0; i < sendCount[ProcRank]; i++)
    {
      residual = fLocal[i];
      for (j = 0; j < n; j++)
      {
        residual -= ALocal[i][j] * x[j];
      }
      residual /= ALocal[i][i + displs[ProcRank]];
      // вычислям часть бесконечной нормы невязки на кажддом процессе
      LocalResidualNorm = std::max(LocalResidualNorm, fabs(residual));
      // вычисляем новое приближение
      x[i + displs[ProcRank]] += residual;
    }
    // вычисляем бесконечную норму из частей, посчитанных на каждос процессе
		MPI_Allreduce(&LocalResidualNorm, &ResidualNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    // союираем части вектора x, посчитанные на каждом процессе и результат прописываем в вектор x на куаждом процессе
    MPI_Allgatherv(MPI_IN_PLACE, 0, 0, x, sendCount, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    // проверяем условие выхода
    if (ResidualNorm < eps)
    {
      // завершаем цикл while(1)
      break;
    }
  }

	if (ProcRank == 0)
  {
    printf("\nSolution: \n");
    for (i = 0; i < n; i++)
    {
      printf("x[%2d] = %lf\n", i, x[i]);
    }
	}

  // освобождаем выделенную динамическую память
  if (ProcRank == 0)
  {
    free_matrix(A);
    free(f);
  }
  free_matrix(ALocal);
  free(fLocal);
  free(x);
  free(sendCount);
  free(displs);

  MPI_Finalize();

	return(0);
}