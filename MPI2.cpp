// MPI_Lab1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
using namespace std;

#define BYTE unsigned char
#define root_proc 0

void show_matrix(BYTE* matrix, const int n, const int m)
{
	if((n<10)&(m<10))
	{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			printf("%5d ", matrix[i * m + j]);
		printf("\n");
	}
	printf("\n");
	}
}

void min_max_search(const BYTE *vector, const int n, BYTE &minimum, BYTE &maximum)
{
	BYTE min = 255, max = 0;
		for (int i = 0; i < n; ++i)
        {
                if (vector[i] < min) min = vector[i];
                if (vector[i] > max) max = vector[i];
        }
        minimum = min;
        maximum = max;
}

int main(int argc, char* argv[])
{
		BYTE* matrix = 0;
		BYTE* matrixp = 0;
		BYTE* matrixl = 0;
		int n, m;
		
		BYTE min,max;
		BYTE max_min_arr[2];
		
		int proc_num, proc_rank;
		int *send_counts, *displs;
		BYTE *recieve_buffer;

		int chunk_size, rem;
		double parallel_start, parallel_end;
		double linear_start, linear_end;
		MPI_Init(&argc, &argv);

		MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
		
		if (proc_rank == root_proc)
		{
			cin>>n;
			cin>>m;
			matrix = new BYTE[n * m];
 			srand(MPI_Wtime() * 1000);
			for (int i = 0; i < n * m; ++i)
				matrix[i] = rand() % 255;
			show_matrix(matrix, n, m);
			matrixl = matrix;
			matrixp = matrix;
			parallel_start = MPI_Wtime();
		}
		
		
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

		send_counts = new int[proc_num];
		displs = new int[proc_num];
		chunk_size = n * m / proc_num;
		for (int i = 0; i < proc_num; ++i)
		{
			send_counts[i] = chunk_size;
			displs[i] = i * chunk_size;
		}
		send_counts[proc_num - 1] += (n * m) % proc_num; 

		recieve_buffer = new BYTE[send_counts[proc_rank]];
		MPI_Scatterv(matrixp, send_counts, displs, MPI_BYTE, recieve_buffer, send_counts[proc_rank], MPI_BYTE, root_proc, MPI_COMM_WORLD);
		
		min_max_search(recieve_buffer, send_counts[proc_rank],max_min_arr[0],max_min_arr[1]);

		BYTE *temp_buffer = new BYTE[2 * proc_num];

		MPI_Allgather(max_min_arr, 2, MPI_BYTE, temp_buffer, 2, MPI_BYTE, MPI_COMM_WORLD);

		min_max_search(temp_buffer, 2 * proc_num, max_min_arr[0], max_min_arr[1]);
		delete temp_buffer;
		
		for (int i = 0; i < send_counts[proc_rank]; ++i)
			recieve_buffer[i] = (recieve_buffer[i] - max_min_arr[0]) * 128 / (max_min_arr[1] - max_min_arr[0])+64;
		

		MPI_Gatherv(recieve_buffer, send_counts[proc_rank], MPI_BYTE,matrixp, send_counts, displs, MPI_BYTE, root_proc, MPI_COMM_WORLD);

		if (proc_rank == root_proc)
		{
			parallel_end = MPI_Wtime();
			show_matrix(matrixp, n, m);
			printf("parallel version\n");
			printf("%f time spent \n\n",parallel_end - parallel_start);

			linear_start = MPI_Wtime();
			BYTE max_l,min_l;
			min_max_search(matrixl,n*m, min_l, max_l);
			for (int i = 0; i < n*m; ++i)
				matrixl[i] = (matrixl[i] - min_l) * 128 / (max_l - min_l)+64;
			show_matrix(matrixl, n, m);
			linear_end = MPI_Wtime();

			printf("linear version\n");
			printf(" %f time spent\n", linear_end - linear_start);

			printf("time boost is %f\n", linear_end - linear_start - parallel_end + parallel_start);
		}

		MPI_Finalize();
	return 0;
}


