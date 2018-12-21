#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <mpi.h>
#include <cmath>

#include "equation.h"


int main(int argc, char ** argv)
{
	int errCode = 0;

	if ((errCode = MPI_Init(&argc, &argv)) != 0)
	{
		std::cout << "MPI initialisation error" << std::endl;
		return errCode;
	}

	// Переменные для MPI:  общее число процессов и номер текущего процесса
    int size, rank;
	// Число уравнений
	int N;
	// Точность
	double eps;

	// Система уравнений
	equation* e = new equation;
	
	// Число процессов
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	// Номер текущего процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		

	// Ветвь менеджера
	if (rank == 0)
	{
		std::cout << "Input amount of the equations N = ";
		std::cin >> N;

		std::cout << "Input order of precision: eps = 10^-";
		std::cin >> eps;
		eps = pow(10., -eps);

		e->eqGenerate(N, eps);		
	}

	// Передаем параметры задачи от 0-го процесса текущему
	// Размерность системы
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Точность вычислений
	MPI_Bcast(&e->EPS, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Максимальное число итераций
	MPI_Bcast(&e->MAXIT, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Макисмальная невязка
	MPI_Bcast(&e->MAXRES, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Создаем массивы для распределения данных
	int *sendcountsMatrix, *sendcountsVector; // Количество 
	int *displsMatrix, *displsVector; // Сдвиг

	sendcountsMatrix = new int[size];
	sendcountsVector = new int[size];
	displsMatrix = new int[size]; 
	displsVector = new int[size];

	if (rank == 0)
	{
		for (int rank_i = 0; rank_i < size; rank_i++)
		{
			// Заполняем массивы с количеством передаваемых данных: каждому процессу передается sendcountsVector[i] строк матрицы и N*sendcountsVector[i] элементов матрицы
			sendcountsVector[rank_i] = N / size;
			if (rank_i < N - size * sendcountsVector[rank_i])
			{
				sendcountsVector[rank_i]++;
			}
			sendcountsMatrix[rank_i] = N * sendcountsVector[rank_i];

			// Заполняем массив сдвигов: сдвиг текущего = сдвиг предыдущего + число передаваемых предыдущему элементов
			if (rank_i == 0)
			{
				displsVector[rank_i] = 0;
				displsMatrix[rank_i] = 0;
			}
			else
			{
				displsVector[rank_i] = displsVector[rank_i - 1] + sendcountsVector[rank_i - 1];
				displsMatrix[rank_i] = displsMatrix[rank_i - 1] + sendcountsMatrix[rank_i - 1];
			}
		}
	}
	

	MPI_Bcast(sendcountsMatrix, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(sendcountsVector, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displsMatrix, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displsVector, size, MPI_INT, 0, MPI_COMM_WORLD);


	// Текущее приближение
	double *X = new double[N];
	for (int i = 0; i < N; i++)
	{
		X[i] = 0.;
	}

	// Рассылка нужных частей (по sendcounts и displs) матрицы A и вектора R
	double * A_ = new double[sendcountsMatrix[rank]];
	double * R_ = new double[sendcountsVector[rank]];
	MPI_Scatterv(e->A, sendcountsMatrix, displsMatrix, MPI_DOUBLE, A_, sendcountsMatrix[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(e->R, sendcountsVector, displsVector, MPI_DOUBLE, R_, sendcountsVector[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// Точность решения на данном шаге
	double precis = e->MAXRES;
	// Номер итерации
	int iter = 0;

	while ((precis > e->EPS) && (iter <= e->MAXIT))
	{
		// Невязка в рамках данного процесса
		double precis_ = -1;
		// Невязка для i-го элемента
		double Ri_;
	
		// Находим часть нового приближения
		for (int i = 0; i < sendcountsVector[rank]; i++)
		{
			// Вычисляем невязку предыдущей аппроксимации в новом уравнении
			Ri_ = R_[i];
			for (int j = 0; j < N; j++)
			{
				Ri_ -= A_[i * N + j] * X[j];
			}
			Ri_ = Ri_ / A_[i*N + i + displsVector[rank]];
			
			// Вычисляем новое приближение
			X[displsVector[rank] + i] += Ri_;
			// Вычисляем точность для процесса с рангом rank
			precis_ = fmax(fabs(Ri_), precis_);

		}
		// Составляем общий X из посчитанных частей
		MPI_Allgatherv(MPI_IN_PLACE, 0, 0, X, sendcountsVector, displsVector, MPI_DOUBLE, MPI_COMM_WORLD);

		// Вычисляем точность для всех процессов на данной итерации
		MPI_Allreduce(&precis_, &precis, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		// Проверяем, не расходится ли решение
		if (precis > e->MAXRES)
		{
			break;
		}

		// Новая итерация
		iter = iter + 1;
	}

	
	// Ветвь менеджера
	if (rank == 0)
	{

		if (precis < e->EPS)
		{ 
			// Решение найдено
			std::cout << "The solution is found after " << iter << " iterations:" << std::endl;
			printVector(X, N, "cmd");
		}
		else if (precis < e->MAXRES)
		{
			// Превышено максимальное число итераций, но решение не найдено с заданной точностью
			std::cout << "Max iterations limit (" << e->MAXIT << ") is exceeded. Current residual is " << precis << " for x ="<< std::endl;
			printVector(X, N, "cmd");
		}
		else
		{
			// Метод разошелся
			std::cout << "Local residue exceeded max value (" << e->MAXRES <<") after "<< iter <<" iterations." << std::endl;
		}	
		
		delete e;
	}
	
	delete R_;
	delete A_;

	delete sendcountsVector;
	delete sendcountsMatrix;
	delete displsVector;
	delete displsMatrix;

	delete X;

    MPI_Finalize();

	return 0;
}

