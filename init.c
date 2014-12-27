#include <omp.h>

#define N 100000

int main(int argc, char const *argv[])
{

int a[N];
int i;


for (i=0; i<N; i++)
	a[i] = i;


#pragma omp parallel
{
	int sum, i;
	int tid = omp_get_thread_num();
	for (i=tid*N/8; i<(tid+1)*N/8; i++)
		sum += a[i];
}
	return 0;
}
