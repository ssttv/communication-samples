#include <omp.h>

#define N 1024

int main(int argc, char const *argv[])
{
	int i;
	int a[N];
	int sum = 0;

	#pragma omp parallel for reduction(+:sum)
	for (i=0; i < N; i++) {
		sum += a[i];
	}

	return 0;
}


