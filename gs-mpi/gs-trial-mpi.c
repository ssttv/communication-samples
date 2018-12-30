#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Метод Гаусса-Зейделя на основе mpi.h 
// Команда для компиляции: mpicc -o gs-trial-mpi.out gs-trial-mpi.c -lm


// Я просил в этой задаче сделать ленточное распределение данных СЛАУ по процессам. 
// Исходную матрицу и правую часть СЛАУ нужно было задать только на нулевом процессе, а затем раздать каждому процессу 
// с помощью функции MPI_Scatterv часть строк матрицы и правой части СЛАУ. То есть на каждом процессе метод Зейделя 
// должен использовать только свою часть СЛАУ. При этом в методе Зейделя распараллеливать нужно цикл for (i = 0; i <= n; i++),
// а не вложенный в него цикл по элементам строки. То есть межпроцессорный обмен нужно производить только один раз в конце 
// каждой итерации метода Зейделя. У Вас же межпроцессорные обмены на каждой итерации стоят в цикле for (i = 0; i <= n; i++), 
// что очень затратно. При этом для межпроцессорного обмена в методе Зейделя удобно использовать MPI_Allgatherv для сборки нового 
// вектора приближения и MPI_Allreduce для вычисления нормы невязки для условия выхода.


double *x;
int ProcNum; 
int ProcRank; 
void Process(int, double);
int main(int agrc, char* argv[])
{
    int n, i;
	double eps;
	FILE *Out;

	MPI_Init(&agrc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		printf("n = ");
		if (scanf("%d", &n) > 0) printf("OK\n");
		printf("eps = ");
		if (scanf("%lf", &eps) > 0) printf("OK\n");
		printf("\n");
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	x = (double*)malloc((n+1)*sizeof(double));
	for (i = 0; i <= n; i++)
		x[i] = 0;

	Process(n, eps);

	//Запись в файл
	if (ProcRank == 0) {
		Out = fopen("output.txt", "w");
		for (i = 0; i <= n; ++i) {
			fprintf(Out, "%lf\n", x[i]);
		}
		fclose(Out);
	
		printf("\nResult written in file output.txt\n");
	}

	//Завершение MPI
	MPI_Finalize();

	return(0);
}

void Process(int n, double eps)
{
    int i, j, counter = 0;
    double *a, *b, *c, *f, *p;
    double *prev_x;
    double **A, **C;
    double norm;
    double ProcSum, TotalSum;

    a = (double*)malloc((n+1)*sizeof(double));
    b = (double*)malloc((n+1)*sizeof(double));
    c = (double*)malloc((n+1)*sizeof(double));
    f = (double*)malloc((n+1)*sizeof(double));
    p = (double*)malloc((n+1)*sizeof(double));

    prev_x = (double*)malloc((n+1)*sizeof(double));

    A = (double**)malloc((n+1)*sizeof(double*));
    for (i = 0; i <= n; ++i)
        A[i] = (double*)malloc((n+1)*sizeof(double));

    C = (double**)malloc((n+1)*sizeof(double*));
    for (i = 0; i <= n; ++i)
        C[i] = (double*)malloc((n+1)*sizeof(double));

    //Заполняем столбцы данных
    b[0] = 1.;
    c[0] = 0.;
    f[0] = 1.;
    p[0] = 1.;

    for (i = 1; i < n; i++) {
        a[i] = 1.;
        b[i] = -2.;
        c[i] = 1.;
        f[i] = 2./(i*i + 1);
        p[i] = 2.;
    }

    f[n] = -n/3.;
    p[n] = 1.;

    //Формируем матрицу коэффициентов
    for (i = 0; i <= n; i++)
    for (j = 0; j <= n; j++)
        A[i][j] = 0.;

    A[0][0] = b[0]; A[0][1] = c[0];

    for (i = 1; i < n; i++) {
        A[i][i] = b[i];
        A[i][i+1] = c[i];
        A[i][i-1] = a[i];
    }

    for (j = 0; j <= n; j++)
        A[n][j] = p[j];

    //Формируем расчетную матрицу
    for (i = 0; i <= n; i++)
    for (j = 0; j <= n; j++) {
        if (i == j) {
            C[i][j] = 0;
        }
        else {
            C[i][j] = -A[i][j]/A[i][i];
        }
    }

	//На каждом процессе будет вычисляться частичная сумма длиной k
	int k = (n+1) / ProcNum;
	//То есть будут суммироваться только те элементы, которые лежат между i1 и i2
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = n+1;

    if (ProcRank == 0) printf("Calculating started\n\n");

    //Взять бесконечную норму погрешности и вычислить невязку на основном процесе
    while(1) {
        for (i = 0; i <= n; i++)-+++++            prev_x[i] = x[i];

        for (i = 0; i <= n; i++) {
            ProcSum = 0.0;

            for (j = 0; j < i; j++)
                if ((i1 <= j) && (j < i2))
			ProcSum += C[i][j]*x[j];

            for (j = i; j <= n; j++)
                if ((i1 <= j) && (j < i2))
			ProcSum += C[i][j]*prev_x[j];

		TotalSum = 0.0;
		//Сборка частичных сумм ProcSum на процессе с рангом 0
		MPI_Reduce(&ProcSum, &TotalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		x[i] = TotalSum + f[i]/A[i][i];
        }

        //Считаем невязку
	if (ProcRank == 0) {
        	norm = 0;
        	for (i = 0; i <= n; i++)
        	    norm += (x[i] - prev_x[i])*(x[i] - prev_x[i]);
        	norm = sqrt(norm);
	
        	counter++;
	
        	printf("%2d.   %lf\n", counter, norm);
	}
	MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (norm < eps) return;
    }
}