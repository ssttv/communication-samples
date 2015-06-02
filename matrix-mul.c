#include <omp.h>

#define N 1024

int A[N][N], B[N][N], C[N][N];

int main(int argc, char const *argv[])
{
	int i=0, j=0, k=0;

	#pragma omp parallel for
	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			for(k=0; k<N; k++) {
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}

	return 0;
}


