
#include "stdafx.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <random>
#include <ctime>
#include <cstdlib>
#include <conio.h>
using namespace std;

int** create_matrix(int col, int row) {
	int** matrix = 0;	
	matrix = new int*[row];
	for (int i = 0; i < row; i++)
		matrix[i] = new int[col];

	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			matrix[i][j] = (int)(rand() % 1000) - 500;
	return matrix;
}

int* convert_matrix(int** matrix, int col, int row) {
	if (matrix) {
		int* vector = new int[col * row];

		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				vector[i * col + j] = matrix[i][j];
		return vector;
	}
	else return NULL;
}

void show_matrix(int** matrix, int col, int row) {
	if (matrix) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++)
				cout << matrix[i][j] << " ";
			cout << endl;
		}
		cout << endl;
	}
}

void show_matrix_vector(int* vector, int col, int row) {
	if (vector) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++)
				cout << vector[i * col + j] << " ";
			cout << endl;
		}
		cout << endl;
	}
}

void show_vector(int* vector, int count) {
		for (int i = 0; i < count; i++)
			cout << vector[i] << " ";
		cout << endl;
}

int findMax(int* vector, int count)
{
	int max = 0;
		max = vector[0];
		for (int i = 0; i < count; ++i)
			if (max < vector[i])
				max = vector[i];
	return max;
}

int main(int argc, char* argv[])
{	
	srand((unsigned)time(NULL));
	// Matrix
	int m_col = 10, m_row = 2;
	int** matrix_matrix = 0;
	int* matrix_vector = 0;
	// Maximum
	int* max_arr = 0;
	int** max_in_row = 0;	
	// Time
	double tline_start, tline_end, tparal_start, tparal_end;
	// MPI
	int proc_rank, proc_size;
	int *send_counts, *displs, *recieve_buffer, *gather_buffer;
	int chunk, rem;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
	// Root Proc gen matrix
	if (proc_rank == 0){
		cout << " row count: "; cin >> m_row;
		cout << " column count: "; cin >> m_col;
		cout << " " << endl;		
		matrix_matrix = create_matrix(m_col, m_row);
		matrix_vector = convert_matrix(matrix_matrix, m_col, m_row);
		max_in_row = new int*[m_row];
		for (int i = 0; i < m_row; i++)
			max_in_row[i] = new int[proc_size];
		max_arr = new int[m_row];
		tparal_start = MPI_Wtime();
	}

	MPI_Bcast(&m_col, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m_row, 1, MPI_INT, 0, MPI_COMM_WORLD);

	send_counts = new int[proc_size];
	displs = new int[proc_size];
	chunk = m_col / proc_size;
	rem = m_col % proc_size;
	for (int i = 0; i < rem; i++)
	{
		send_counts[i] = chunk + 1;
		displs[i] = i * (chunk + 1);
	}
	for (int i = rem; i < proc_size; i++)
	{
		send_counts[i] = chunk;
		displs[i] = rem + i * chunk;
	}
	
	recieve_buffer = new int[send_counts[proc_rank]];	
	gather_buffer = new int[proc_size];
	int maximum, k = 0;

	for (k = 0; k < m_row; k++) 
	{		
		MPI_Scatterv(matrix_vector, send_counts, displs, MPI_INT, recieve_buffer, send_counts[proc_rank], MPI_INT, 0, MPI_COMM_WORLD);
		
		for (int i = 0; i < proc_size; i++)
			displs[i] += m_col;

		maximum = findMax(recieve_buffer, send_counts[proc_rank]);			

		MPI_Gather(&maximum, 1, MPI_INT, gather_buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);		
		
		if (proc_rank == 0)
		{
			for (int i = 0; i < proc_size; i++)
				max_in_row[k][i] = gather_buffer[i];
		}
	}

	if (proc_rank == 0)
	{
		for (int i = 0; i < m_row; i++)
			max_arr[i] = findMax(max_in_row[i], proc_size);

		tparal_end = MPI_Wtime();

		cout <<endl << "  PARALLEL VERSION  " << endl;
		if (m_row < 20) {
			for (int i = 0; i < m_row; i++)
				cout << " maximum value in row " << i << " is " << max_arr[i] << endl;
		}
		cout << " Parallel time: " << tparal_end - tparal_start << " ms" << endl;
		
		tline_start = MPI_Wtime();
		int *linear_max = new int[m_col];
		for (int i = 0; i < m_row; i++)			
			linear_max[i] = findMax(matrix_matrix[i], m_col);
		tline_end = MPI_Wtime();
		
		cout << endl << "  LINEAR VERSION  " << endl;
		if (m_row < 20) {
			for (int i = 0; i < m_row; i++)
				cout << " maximum value in row " << i << " is " << linear_max[i] << endl;
		}
		cout << " Linear time: " << tline_end - tline_start << " ms" << endl;

		cout << " Boost: " << ((tline_end - tline_start) - (tparal_end - tparal_start)) << endl;
	}
	
	MPI_Finalize();		
    return 0;
}

