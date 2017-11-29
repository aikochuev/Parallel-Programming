#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <math.h>
#include <string.h>

using namespace std;

#define ROOT 0

int procRank, procCount;
int *procPivotIter;

void MultiMatrix(double *A, double *B, double *C, const int n)  	
{    	
    for(int i = 0; i < n; ++i)      	
        for(int j = 0; j < n; ++j)        	
            for(int k = 0; k < n; ++k)   	
                C[i * n + j] += A[i * n + k] * B[k * n + j];    	
}

struct PivotChoice {
	double max;
	int procRank;
};

void swap(double *v1, double *v2) {
	double tmp = *v1;
	*v1 = *v2;
	*v2 = tmp;
}

void swapRows(double *r1, double *r2, const int n) {
	int i;

	for (i = 0; i < n; ++i)
		swap(&r1[i], &r2[i]);
}

double *createVec(const int n) {
	double *vec;
	int i;

	vec = new double[n];

	srand(time(0));

	for (i = 0; i < n; ++i)
		vec[i] = 1 + rand() % 10;
	
	return vec;
}

double *createVecElem(const int n) {
	double *vec;
	int i;
	int h=0;
	vec = new double [n * n];

	for (i = 0; i < n*n; ++i)
	{
		if(i==h*(n+1))
		{
			vec[i]=1;
			h++;
		}
		else
			vec[i]=0;
	}
	return vec;
}

void printMatrix(double *a, double *b, const int rowsCount, const int colsCount) {
	int i, j,h;

	for (i = 0; i < rowsCount; ++i) {
		for (j = 0; j < colsCount; ++j)
			printf("%3.2f ", a[i * colsCount + j]);
		for (h=0;h<colsCount;++h)
			printf("%3.2f ", b[i * colsCount + h]);
		printf("\n");
	}
}

void divideVec(double *vec, const int n, const double div) 
{
	int i;

	for (i = 0; i < n; ++i)
		vec[i] /= div;
}

void subVecs(double *minuend, double *sub, const int n, const double mult) 
{
	int i;

	for (i = 0; i < n; ++i)
		minuend[i] -= (mult * sub[i]);
}

void cpyVec(double *source, double *dest, const int n) 
{
	int i;

	for (i = 0; i < n; ++i)
		dest[i] = source[i];
}

void LineAlgoritm(double *MatrixA,double *MatrixB,int n)
{
	double max;
	int pivotPos; 
	double *pivotRowA = new double [n]; 
	double *pivotRowB = new double [n]; 

	procPivotIter = new int [n];
	for (int i = 0; i < n; ++i)
		procPivotIter[i] = -1;

	for (int i = 0; i < n; i++)
	{
		 max = 0.0;
		 for (int j = 0; j < n; j++)
			if (max < fabs(MatrixA[j * n + i]) && procPivotIter[j] == -1) 
			{
				max = fabs(MatrixA[j * n + i]);
				pivotPos = j;
			}
		//cout<<"max"<<max<<endl;
		 if (max != 0.0) 
		 {
			double divA	= MatrixA[pivotPos * n + i];        //То делем соответствующую строку на ведущий элемент
			//cout<<"divA"<<divA<<endl;
			procPivotIter[pivotPos] = i;    //Запоминаем, на какой итерации внешнего цикла строка была использована
			divideVec(&MatrixA[pivotPos * n], n, divA);     //Делим на ведущий элемент
			divideVec(&MatrixB[pivotPos * n], n, divA);
			cpyVec(&MatrixA[pivotPos * n], pivotRowA, n);
			cpyVec(&MatrixB[pivotPos * n], pivotRowB, n);//Копируем
		 }
		 for (int j = 0; j < n; ++j)//Эта функция вычитает отправленную строку
			if (procPivotIter[j] == -1) 
			{
				double multA = MatrixA[j * n + i];
				/*double multB = MatrixB[j * n + i];
				cout<<"multA"<<multA<<endl;
				cout<<"multB"<<multB<<endl;*/

				subVecs(&MatrixA[j * n], pivotRowA, n, multA);
				subVecs(&MatrixB[j * n], pivotRowB, n, multA);
			}
		/*cout<<"Output matrix Line: "<<endl;
				printMatrix(MatrixA, MatrixB, n, n);*/
	}

	double *pivotElemA = new double [n];
	double *pivotElemB = new double [n];

	for (int i = n - 1; i >= 0; --i) 
	{
		pivotPos = 0;
		for (int j = 1; j < n; ++j)
			if (procPivotIter[pivotPos] < procPivotIter[j])
				pivotPos = j;

		max = procPivotIter[pivotPos];
		//cout<<"max"<<max<<endl;
		procPivotIter[pivotPos] = -1;
		cpyVec(&MatrixB[pivotPos * n], pivotElemB, n);
		cpyVec(&MatrixA[pivotPos * n], pivotElemA, n);

		for (int j = 0; j < n; ++j)
			if (procPivotIter[j] != -1) 
			{
				double multA = MatrixA[j * n + i];
				//cout<<"multA"<<multA<<endl;
		
				subVecs(&MatrixB[j * n], pivotElemB, n, multA);
				subVecs(&MatrixA[j * n], pivotElemA, n, multA);
			}
			/*cout<<"Output matrix Line: "<<endl;
				printMatrix(MatrixA, MatrixB, n, n);*/
	}

	int k = 0;
	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)
			if (MatrixA[i * n + j] == 1.0) 
			{
				swapRows(&MatrixA[i * n], &MatrixA[k * n], n);
				swapRows(&MatrixB[i * n], &MatrixB[k * n], n);
				++k;
			}
}

void parallelDirectFlow(double *procRowsA, double *procRowsB, const int n, const int rowsCount) {
        double max;             //Максимальное значение среди строк, которые выделили этому процессу
        int pivotPos;   //Соответственно номер этой строки
        struct PivotChoice procPivot, pivot;    //procPivot нужна, чтобы обменяться с другими процессами и узнать максимум, в pivot будет храниться этот максимум
        double *pivotRowA = new double [n];  //Здесь будет храниться целая строка, которую будем вычитать из других
		double *pivotRowB = new double [n]; 
 
        for (int i = 0; i < n; ++i) 
		{    //Основной цикл, нужно n шагов, чтобы алгоритм Гаусса-Жордана работал.
                max = 0.0;      //Находим максимум
                for (int j = 0; j < rowsCount; ++j)
                        if (max < fabs(procRowsA[j * n + i]) && procPivotIter[j] == -1) 
						{
                                max = fabs(procRowsA[j * n + i]);
                                pivotPos = j;
                        }
                       
                procPivot.max = max;    //Присваиваем найденный максимум переменной, которой обменяемся с другими процессами
                procPivot.procRank = procRank;
 
                MPI_Allreduce(&procPivot, &pivot, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);       //Выполняет операцию поиска максимум среди procPivot всех процессов и ложит его в pivot
 
                if (pivot.max != 0.0) 
				{
                        if (procRank == pivot.procRank) 
						{    //Если максимум достигнут для этого процесса
                                double divA	= procRowsA[pivotPos * n + i];        //То делем соответствующую строку на ведущий элемент

                                procPivotIter[pivotPos] = i;    //Запоминаем, на какой итерации внешнего цикла строка была использована
 
                                divideVec(&procRowsA[pivotPos * n], n, divA);     //Делим на ведущий элемент
								divideVec(&procRowsB[pivotPos * n], n, divA);
 
                                cpyVec(&procRowsA[pivotPos * n], pivotRowA, n);
								cpyVec(&procRowsB[pivotPos * n], pivotRowB, n);//Копируем
                        }
 
                        MPI_Bcast(pivotRowA, n , MPI_DOUBLE, pivot.procRank, MPI_COMM_WORLD); //Отправляем строку, которую надо вычитать другим процессам
						MPI_Bcast(pivotRowB, n , MPI_DOUBLE, pivot.procRank, MPI_COMM_WORLD);
                          
						for (int j = 0; j < rowsCount; ++j)//Эта функция вычитает отправленную строку
							if (procPivotIter[j] == -1) 
							{
								double multA = procRowsA[j * n + i];
								//double multB = procRowsB[j * n + i];

								subVecs(&procRowsA[j * n], pivotRowA, n, multA);
								subVecs(&procRowsB[j * n], pivotRowB, n, multA);
							}
                }
        }
 
        delete [] pivotRowA;
		delete [] pivotRowB;
}

void parallelReverseFlow(double *procRowsA, double *procRowsB, const int n, const int rowsCount) 
{
	int pivotPos;
	struct PivotChoice procPivot, pivot;
	double *pivotElemA = new double [n];
	double *pivotElemB = new double [n];

	for (int i = n - 1; i >= 0; --i) 
	{
		pivotPos = 0;
		for (int j = 1; j < rowsCount; ++j)
			if (procPivotIter[pivotPos] < procPivotIter[j])
				pivotPos = j;

		procPivot.max = procPivotIter[pivotPos];
		procPivot.procRank = procRank;

		MPI_Allreduce(&procPivot, &pivot, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

		if (procRank == pivot.procRank) 
		{
			procPivotIter[pivotPos] = -1;
			cpyVec(&procRowsB[pivotPos * n], pivotElemB, n);
			cpyVec(&procRowsA[pivotPos * n], pivotElemA, n);
		}

		MPI_Bcast(pivotElemA, n, MPI_DOUBLE, pivot.procRank, MPI_COMM_WORLD);
		MPI_Bcast(pivotElemB, n, MPI_DOUBLE, pivot.procRank, MPI_COMM_WORLD);

		for (int j = 0; j < rowsCount; ++j)
			if (procPivotIter[j] != -1) 
			{
				double multA = procRowsA[j * n + i];
				//double multB = procRowsB[j * n + i];
				//procRowsA[j * n + i] = 0.0;
				subVecs(&procRowsB[j * n], pivotElemB, n, multA);
				subVecs(&procRowsA[j * n], pivotElemA, n, multA);
			}
	}
}

int main(int argc, char *argv[]) {
	double *MatrixA = NULL;
	double *MatrixB = NULL;
	double *resMatrixA = NULL;
	double *resMatrixB = NULL;
	double *procRowsA = NULL;
	double *procRowsB = NULL;
	double parallel_start, parallel_end;
	double line_start, line_end;
	int n=atoi(argv[1]);

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);

	if (procRank == ROOT) 
	{
		MatrixA = createVec(n * n);
		MatrixB = createVecElem(n);
		resMatrixA = new double [n * n];
		resMatrixB = new double [n * n];
		parallel_start = MPI_Wtime();
	}
	if (procRank == ROOT && n<8) 
	{
		cout<<"Source matrix: "<<endl;
		printMatrix(MatrixA, MatrixB, n, n);
	}

	int *displs=new int [procCount];
	int *sendCount=new int [procCount];
	int chunk=n/procCount;
	for(int i=0;i<procCount;++i)
	{
		sendCount[i]=chunk*n;
		displs[i]=i*chunk*n;
	}
	sendCount[procCount-1]+=(n%procCount)*n;

	procPivotIter = new int [sendCount[procRank]/n];
	for (int i = 0; i < sendCount[procRank]/n; ++i)
		procPivotIter[i] = -1;

	procRowsA=new double [sendCount[procRank]];
	procRowsB=new double [sendCount[procRank]];

	MPI_Scatterv(MatrixA, sendCount, displs, MPI_DOUBLE, procRowsA, sendCount[procRank], MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	MPI_Scatterv(MatrixB, sendCount, displs, MPI_DOUBLE, procRowsB, sendCount[procRank], MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	/*cout<<"Process:"<<procRank<<endl;
	printMatrix(procRowsA, procRowsB, sendCount[procRank]/n, n);*/

	parallelDirectFlow(procRowsA,procRowsB,n, sendCount[procRank]/n);
	//MPI_Barrier(MPI_COMM_WORLD);
	parallelReverseFlow(procRowsA,procRowsB,n,sendCount[procRank]/n);

	
	MPI_Gatherv(procRowsA, sendCount[procRank], MPI_DOUBLE, resMatrixA, sendCount, displs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	MPI_Gatherv(procRowsB, sendCount[procRank], MPI_DOUBLE, resMatrixB, sendCount, displs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	
	if (procRank == ROOT) 
	{
		int k = 0;
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i)
				if (resMatrixA[i * n + j] == 1.0) 
				{
					swapRows(&resMatrixA[i * n], &resMatrixA[k * n], n);
					swapRows(&resMatrixB[i * n], &resMatrixB[k * n], n);
					++k;
				}
		parallel_end = MPI_Wtime();
		double *TEST=new double [n*n]();
		MultiMatrix(MatrixA,resMatrixB,TEST,n);
		if (n < 8) 
		{
		cout<<"Output matrix: "<<endl;
				printMatrix(resMatrixA, resMatrixB, n, n);
		}
		printf("%f time spent PARALLEL \n\n",parallel_end - parallel_start);
		line_start = MPI_Wtime();
		LineAlgoritm(MatrixA,MatrixB,n);
		line_end = MPI_Wtime();
		double Delta=0.0;
		double maxDelta=0.0;
		for(int i=0;i<n*n;i++)
		{
			Delta=resMatrixB[i]-MatrixB[i];
			if(fabs(Delta)>fabs(maxDelta))
				maxDelta=Delta;
		}
		printf("%g \n", maxDelta);
		if (n < 8) 
		{
		cout<<"Output matrix Line: "<<endl;
				printMatrix(MatrixA, MatrixB, n, n);
		cout<<"TEST matrix: "<<endl;
				printMatrix(TEST, TEST, n, n);
		}
		printf("%f time spent LINE \n\n",line_end - line_start);
		delete [] sendCount;
		delete [] displs;
	}
	MPI_Finalize();
}