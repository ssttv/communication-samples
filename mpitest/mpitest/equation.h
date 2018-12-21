#pragma once
#include <string>
using std::string;

class equation
{
public:
	// Число уравнений
	int N;
	// Точность вычислений
	double EPS;
	// Максимальное число итераций
	int MAXIT = 200;
	// Максимальная невязка
	double MAXRES = 1e10;

	// Матрица системы
	double *A = NULL;
	// Столбец правой части
	double *R = NULL;

	// Сгенерировать матрицу и вектор правой части системы
	void eqGenerate(int n, double eps);
};

// Напечатать произвольный вектор в файл fname
int printVector(double *v, int N, string fname);
