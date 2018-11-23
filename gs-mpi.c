#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
 
/*
 * Команда для компиляции: mpicc -o gs-mpi gs-mpi.c
 */

void test_system
 ( int n, double **A, double *b );

void run_gauss_seidel_method
 ( int p, int n, double **A, double *b,
   double epsilon, int maxit,
   int *numit, double *x );

/*
 * Тут запускается метод Гаусса-Зейдееля для A*x = b.
 *
 * Вход:
 *   p        число потоков;
 *   n        размерность матрицы;
 *   A        матрица n на n где A[i][i] /= 0 (диагональная);
 *   b        n-размерный вектор;
 *   epsilon  мера погрешность;
 *   maxit    макисмальное число итераций;
 *   x        начальный вектор для отдельной итерации.
 *
 * ВЫХОД :
 *   numit    число пройденных итераций;
 *   x        аппроксимация решения A*x = b. 
 * */



int main ( int argc, char *argv[] )
{
   int n,p,i;
   if(argc > 1)
   {
      n = atoi(argv[1]);
      p = (argc > 2) ? atoi(argv[2]) : 1;
   }
   else
   {
      printf("Введите число элементов в строках и столбцах матрицы : ");
      scanf("%d",&n);
      printf("Введите число потоков : ");
      scanf("%d",&p);
   }
   omp_set_num_threads(p);
   {
      double *b;
      b = (double*) calloc(n,sizeof(double));
      double **A;
      A = (double**) calloc(n,sizeof(double*));
      for(i=0; i<n; i++)
         A[i] = (double*) calloc(n,sizeof(double));
      test_system(n,A,b);
      double *x;
      x = (double*) calloc(n,sizeof(double));
      for(i=0; i<n; i++) x[i] = 0.0;
      double eps = 1.0e-4;
      int maxit = 2*n*n;
      int cnt = 0;
      run_gauss_seidel_method(p,n,A,b,eps,maxit,&cnt,x);
      printf("Вычислены %d итераций\n",cnt);
      double sum = 0.0;
      for(i=0; i<n; i++)
      {
         double d = x[i] - 1.0;
         sum += (d >= 0.0) ? d : -d;
      }
      printf("Погрешность : %.3e\n",sum);
   }
   return 0;
}

void test_system
 ( int n, double **A, double *b )
{
   int i,j;
   for(i=0; i<n; i++)
   {
      b[i] = 2.0*n;
      for(j=0; j<n; j++) A[i][j] = 1.0;
      A[i][i] = n + 10.0;
   }
}

void run_gauss_seidel_method
 ( int p, int n, double **A, double *b,
   double epsilon, int maxit,
   int *numit, double *x )
{
   double *dx;
   dx = (double*) calloc(n,sizeof(double));
   int i,j,k,id,jstart,jstop;

   int dnp = n/p;
   double dxi; 
   
   double sum = 0.0;
   #pragma omp parallel \
    shared(A, x, sum) private(k, dxi)

   for(k=0; k<maxit; k++)
   {
      #pragma omp barrier
      sum = 0.0;
      double localSum = 0.0;
      #pragma omp for   
      for(i=0; i<n; i++)
      { 
        localSum = 0.0;
         dxi = b[i];
         for(j=0; j<n; j++)
            dxi -= A[i][j]*x[j];
         dxi /= A[i][i];
         x[i] += dxi;
         localSum += ( (dxi >= 0.0) ? dxi : -dxi);
      }
      #pragma omp critical
      sum += localSum;

      #pragma omp barrier

      printf("%4d : %.3e\n",k,sum);
      if(sum <= epsilon) break;
   }
   //*numit = k+1;
   free(dx);
}
