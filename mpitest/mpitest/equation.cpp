#include "equation.h"

#include <string>
#include <iostream>
#include <fstream>

#include <ctime>
#include <cmath>


using std::ofstream;
using std::string;
using std::cout;

// ������������� ������� � ������ ������ ����� ������� (��������� A � R)
void equation::eqGenerate(int n, double eps)
{
	N = n;
	EPS = eps;
	A = new double[N*N];
	R = new double[N];

	// ������� ������� ������������ ��� A = Q*B*Q', 
	// ��� B - ��������� ������������ �������, 
	// Q - ������������� �������, ���������� ��������������� ����������� �� ���������� ������� v: Q = I-vv'/(v'v)
	double *Q = new double[N*N];
	double *B = new double[N];
	double *v = new double[N];
	double v_norm = 0;
	// ������ ������� ������������ �� ������������ ������� A � ������� ������� x0 = 1..N: b = A*x0
	double *x0 = new double[N];

	//srand(time(NULL));
	srand(1);

	for (int i = 0; i < N; i++)
	{
		// ���������� v
		v[i] = double(1 + rand() % 10);
		v_norm += v[i] * v[i];
		// ���������� B
		B[i] = double(1 + rand() % 10);
		// ���������� x0
		x0[i] = double(i + 1);
		// �������������� R
		R[i] = 0.;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// ���������� Q
			Q[i*N + j] = - 2.* v[i] * v[j]/v_norm;
			// �������������� A
			A[i*N + j] = 0.;
			
		}
		Q[i*N + i] += 1. ;
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				// ��������� A
				A[i*N + j] += Q[i*N+k]*B[k]*Q[j*N+k];
			}
			// ��������� R
			R[i] += A[i*N + j] * x0[j];
		}
	}

	delete[] Q;
	delete[] B;
	delete[] v;
	delete[] x0;
}


// ���������� ������������ ������ � ���� fname
int printVector(double *v, int N, string fname)
{
	if (fname == "cmd")
	{
		for (int i = 0; i < N; i++)
		{
			cout << v[i] << " \n";
		}
		cout << " \n";
		return 0;
	}
	else
	{
		ofstream ofile(fname);

		if (ofile.is_open())
		{
			for (int i = 0; i < N; i++)
			{
				ofile << v[i] << " \n";
			}
			ofile << " \n";
			ofile.close();
			return 0;
		}
		else
		{
			return 1;
		}
	}
}




