#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include "tbb\tbb.h"
#include "tbb\task_scheduler_init.h"
#include <ctime>

using namespace std;
using namespace tbb;
int ThreadNum1;
class MyTask {
	int RowIndex_;
	int ColIndex_;
	int Size_;
	int *pAMatrix_;
	int *pCMatrix_;
	int *pBMatrix_;
	int ThreadNum_;
public:
	MyTask(int ThreadNum, int *pAMatrix, int *pBMatrix, int *pCMatrix, int RowIndex, int ColIndex, int Size) : ThreadNum_(ThreadNum),pAMatrix_(pAMatrix), pBMatrix_(pBMatrix), pCMatrix_(pCMatrix), RowIndex_(RowIndex), ColIndex_(ColIndex), Size_(Size) {}
	void operator()() 
	{
		int GridSize = int(sqrt((int)ThreadNum_));
		int BlockSize = Size_ / GridSize;
		//cout << "task" << name_ << endl;
		for (int iter = 0; iter < GridSize; iter++)
		{
			for (int i = RowIndex_*BlockSize; i < (RowIndex_ + 1)*BlockSize; i++)
				for (int j = ColIndex_*BlockSize; j < (ColIndex_ + 1)*BlockSize; j++)
					for (int k = iter*BlockSize; k < (iter + 1)*BlockSize; k++)
					{
						pCMatrix_[i*Size_ + j] += pAMatrix_[i*Size_ + k] * pBMatrix_[k*Size_ + j];
						//if(ThreadID==0)
						//cout <<"  iter = "<<iter<< "___"<<pAMatrix_[i*Size_+k]<<" * "<< pBMatrix_[k*Size_ + j] <<" = "<<pCMatrix_[i*Size_ + j] << endl;
					}
		}
	}
};
int *GetRandomMatrix(const int n)
{
	int *result = new int[n * n];

	int rand_max = 25;
	for (int i = 0; i < n * n; ++i)
		result[i] = rand() % rand_max;
	return result;
}
void MultiplyMatrices(int *A, int *B, int *C, const int n)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				C[i * n + j] += A[i * n + k] * B[k * n + j];
}
void ParallelResultCalculation(int* pAMatrix, int* pBMatrix, int *pCMatrix, int Size) 
{
	int GridSize = int(sqrt((int)ThreadNum1));
	int BlockSize = Size / GridSize;
	omp_set_num_threads(ThreadNum1);
#pragma omp parallel
	{
		int ThreadID = omp_get_thread_num();
		int RowIndex = ThreadID / GridSize;
		int ColIndex = ThreadID%GridSize;
		for (int iter = 0; iter<GridSize; iter++)
		{
			for (int i = RowIndex*BlockSize; i<(RowIndex + 1)*BlockSize; i++)
				for (int j = ColIndex*BlockSize; j<(ColIndex + 1)*BlockSize; j++)
					for (int k = iter*BlockSize; k < (iter + 1)*BlockSize; k++)
					{
						pCMatrix[i*Size + j] += pAMatrix[i*Size + k] * pBMatrix[k*Size + j];
						/*if(ThreadID==0)
							cout << "ThreadNum" << ThreadID<<"  iter = "<<iter<< "___"<<pAMatrix[i*Size+k]<<" * "<< pBMatrix[k*Size + j] <<" = "<<pCMatrix[i*Size + j] << endl;*/
					}
		}
	}
}
bool Error(int *a, int*b, int n)
{
	for (int i = 0; i < n; i++)
	{
		if (a[i] == b[i])
			return 1;
		else
			return 0;
	}
}
void printMatrix(int *a, int*b, int *c, const int rowsCount, const int colsCount) 
{
	int i, j, h, k;
	if (rowsCount < 25)
	{
		for (i = 0; i < rowsCount; ++i) {
			for (j = 0; j < colsCount; ++j)
				printf("%i ", a[i * colsCount + j]);
			cout << '|';
			for (h = 0; h < colsCount; ++h)
				printf("%i ", b[i * colsCount + h]);
			cout << '|';
			printf("\n");
		}
		for (i = 0; i < rowsCount; ++i) {
			for (k = 0; k < colsCount; ++k)
				printf("%i ", c[i * colsCount + k]);
			cout << '|';
			printf("\n");
		}
	}
}
void main()
{
	srand(time(0));
	cout << omp_get_max_threads() << endl;
	int n;
	cout << "n = ";
	cin >> n;
	cout <<"ThreadNum = ";
	cin >> ThreadNum1;
	int *A = new int[n * n];
	int *B = new int[n * n];
	int *C1 = new int[n * n];
	int *C2 = new int[n * n];
	int *C3 = new int[n * n];
	A = GetRandomMatrix(n);
	B = GetRandomMatrix(n);
	for (int i = 0; i < n * n; ++i)
		C1[i] = 0;
	for (int i = 0; i < n * n; ++i)
		C2[i] = 0;
	for (int i = 0; i < n * n; ++i)
		C3[i] = 0;
	double t_start = 0, t_end = 0, serial = 0, parallel_OpenMP = 0, parallel_TBB = 0;
	//printMatrix(A, B, C1, n, n);
	t_start = omp_get_wtime();
	ParallelResultCalculation(A, B, C1, n);
	t_end = omp_get_wtime();
	parallel_OpenMP = t_end - t_start;
	cout << "Parallel_OpenMP = " << parallel_OpenMP<< endl;
	printMatrix(A, B, C1, n, n);
	//cout << "Acs = " << serial / parallel << endl;
	//printMatrix(A, B, C, n, n);
	t_start = 0;
	t_end = 0;
	task_scheduler_init init(ThreadNum1);
	t_start = omp_get_wtime();
	task_group tg;
		for (int k = 0; k < sqrt((int)ThreadNum1); k++)
			for (int j = 0; j < sqrt((int)ThreadNum1); j++)
				tg.run(MyTask(ThreadNum1, A, B, C3, k, j, n));
	tg.wait();
	t_end = omp_get_wtime();
	parallel_TBB = t_end - t_start;
	cout << "Parallel_TBB = " << parallel_TBB << endl;
	printMatrix(A, B, C3, n, n);
	t_start = 0;
	t_end = 0;
	t_start = omp_get_wtime();
	MultiplyMatrices(A, B, C2, n);
	t_end = omp_get_wtime();
	serial = t_end - t_start;
	printMatrix(A, B, C2, n, n);
	cout << "Serial = " << serial << endl;
	cout << Error(C1, C2, n) << endl;
	cout << Error(C2, C3, n) << endl;
	system("pause");
}

