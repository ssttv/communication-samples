/* L-20 MCS 572 Fri 28 Feb 2014 : gauss_seidel_omp.c
 * With OpenMP we run the method of Gauss-Seidel on a test system A*x = b,
 * where A is diagonally dominant and the exact solution consists
 * of all ones.  The user can provide the dimension and the number
 * of threads at the command line. */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void test_system
 ( int n, double **A, double *b );
/*
 * Given n on entry,
 * On return is an n-by-n matrix A
 * with n+1 on the diagonal and 1 elsewhere.
 * The elements of the right hand side b
 * all equal 2*n, so the exact solution x
 * to A*x = b is a vector of ones. */

void run_gauss_seidel_method
 ( int p, int n, double **A, double *b,
   double epsilon, int maxit,
   int *numit, double *x );
/*
 * Runs the  method of Gauss-Seidel for A*x = b.
 *
 * ON ENTRY :
 *   p        number of threads;
 *   n        the dimension of the system;
 *   A        an n-by-n matrix A[i][i] /= 0;
 *   b        an n-dimensional vector;
 *   epsilon  accuracy requirement;
 *   maxit    maximal number of iterations;
 *   x        start vector for the iteration.
 *
 * ON RETURN :
 *   numit    number of iterations used;
 *   x        approximate solution to A*x = b. */

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
      printf("give the dimension : ");
      scanf("%d",&n);
      printf("Give the number of threads : ");
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
      /* we start at an array of all zeroes */
      for(i=0; i<n; i++) x[i] = 0.0;
      double eps = 1.0e-4;
      int maxit = 2*n*n;
      int cnt = 0;
      run_gauss_seidel_method(p,n,A,b,eps,maxit,&cnt,x);
      printf("computed %d iterations\n",cnt);
      double sum = 0.0;
      for(i=0; i<n; i++) /* compute the error */
      {
         double d = x[i] - 1.0;
         sum += (d >= 0.0) ? d : -d;
      }
      printf("error : %.3e\n",sum);
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
      A[i][i] = n + 1.0;
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

   for(k=0; k<maxit; k++)
   {
      double sum = 0.0;
      for(i=0; i<n; i++)
      {
         dx[i] = b[i];
         #pragma omp parallel \
            shared(A,x) \
            private(id,j,jstart,jstop,dxi)
         {
            id = omp_get_thread_num();
            jstart = id*dnp;
            jstop = jstart + dnp;
            dxi = 0.0;
            for(j=jstart; j<jstop; j++)
               dxi += A[i][j]*x[j];
            #pragma omp critical
               dx[i] -= dxi;
         }
         dx[i] /= A[i][i];
         x[i] += dx[i];
         sum += ( (dx[i] >= 0.0) ? dx[i] : -dx[i]);
      }
      printf("%4d : %.3e\n",k,sum);
      if(sum <= epsilon) break;
   }
   *numit = k+1;
   free(dx);
}
