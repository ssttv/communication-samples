#include <iostream> 
#include <stdlib.h>
#include <cmath>
#include <thread>
#include <pthread.h>
#include <omp.h>
#include <semaphore.h>
#include <stack>
using namespace std;

void Gauss_parallel(int n, int m){
 
    float **matrix = new float *[n]; 
    for (int i = 0; i < n; i++)
        matrix[i] = new float [m];
    // Инициализируем матрицу Гаусса
    // n - число уравнений
    // m - число переменных
    for (int i = 0; i < n; i++)
 
        for (int j = 0; j < m; j++)
        {
            matrix[i][j] = rand()%10;
        }
        // Выводим на экран
        /*cout << "Матрица: " << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
                cout << matrix[i][j] << " ";
            cout << endl;
        }
        cout << endl;*/
 
        float  tmp, xx[10000];
        int i, j, k;
        double dt, timein, timeout;
        // Метод Гаусса, прямой ход
        timein = omp_get_wtime();
        omp_set_num_threads(THREADS); // Число потоков указано в define, их 4
#pragma omp parallel
        {
#pragma omp parallel for firstprivate(j) lastprivate(i) shared(matrix, tmp) 
        for (i = 0; i < n; i++)
        {
            tmp = matrix[i][i];
            for (j = n;j >= i;j--)
                matrix[i][j] /= tmp;
            for (j = i + 1;j < n;j++)
            {
                tmp = matrix[j][i];
                for (k = n;k >= i;k--)
                    matrix[j][k] -= tmp * matrix[i][k];
            }
        }
        //Обратный ход
        xx[n - 1] = matrix[n - 1][n];
#pragma omp parallel for firstprivate(j) lastprivate(i) shared(xx, matrix) 
        for (i = n-2; i >= 0; i--)
        {
            xx[i] = matrix[i][n];
            for (j = i + 1;j < n;j++) 
                xx[i] -= matrix[i][j] * xx[j];
        }
        }
        timeout = omp_get_wtime();
        dt = timeout - timein;
        cout << "Время вычислений (параллельно): " << dt << " секунд" << endl;
 
        for (int i = 0; i < n; i++){
            cout << floor(xx[i] * 1000)/1000. << " ";
            cout << endl;
        }
        delete[] matrix;
}
