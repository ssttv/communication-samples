/*
**  PROGRAM: A simple SPMD producer/consumer program
**
**  PURPOSE: this is just a stupid little program to play around
**  with different ways data is shared between threads.
**
**  HISTORY: Written by Tim Mattson, April 2007.
*/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N        100000


/* function to fill an array with random numbers */
void fill_rand(int tid, int a[])
{
   int i;
   for (i=tid*N/8;i<(tid+1)*N/8;i++)
     a[i] = random();

}

/* function to sum the elements of an array */
int Sum_array(int tid, int a[])
{
   int i;  int sum = 0;
   for (i=tid*N/8;i<(tid+1)*N/8;i++)
    sum += a[i];
   return sum;
}

int main(int argc, char *argv[])
{
  int A[N], sum;
  int flag = 0;

  #pragma omp parallel num_threads(8)
  {
    int tid = omp_get_thread_num();

          if (tid<4) {
           fill_rand(tid, A);
           #pragma omp flush
           flag = 1;
           #pragma omp flush (flag)
        }

         else {
           #pragma omp flush (flag)
           while (flag != 1){
              #pragma omp flush (flag)
           }

           #pragma omp flush
           sum = Sum_array(tid, A);
        }


   }

   printf("the sum is %d \n",sum);

   return 0;
}

