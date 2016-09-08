#include <cstdio>
#include <ctime>

#include <cuda.h>
#include <cuda_runtime.h>


#include "cplxNum.h"
#include "samplingroutine.h"


#ifndef min
#define min(x,y) (x>y?x:y)
#endif


void sampling_routine_test_simpleComplex_float()
{
	const int vector_dim = 4;
	const int sampling_points = 8;

	simpleComplexFloat *host_mat;
	simpleComplexFloat *host_pts;

	int i;

	host_mat = (simpleComplexFloat*)malloc( sizeof( simpleComplexFloat ) * vector_dim * vector_dim );
	host_pts = (simpleComplexFloat*)malloc( sizeof( simpleComplexFloat ) * sampling_points );

	/* definition of test matrix
	m = [ -0.736007-0.486131im   -0.50265+0.303371im   -1.13945-0.492444im    0.521185+1.68789im ;
		   0.058992-0.682839im   -0.24367+0.851119im    1.63444+1.08207im     0.471581-0.283695im ;
           0.427102+1.99722im     0.0900991-2.27228im   0.294868+0.495371im  -1.07317-0.438353im ;
           0.782054-1.08824im     0.261584-0.122245im   0.921061+0.306364im   0.456282+2.01316im ]
	*/

	host_mat[0]  = make_simpleComplexFloat(-0.736007f, -0.486131f);
	host_mat[1]  = make_simpleComplexFloat(-0.50265f,   0.303371f);
	host_mat[2]  = make_simpleComplexFloat(-1.13945f,  -0.492444f);
	host_mat[3]  = make_simpleComplexFloat( 0.521185f,  1.68789f);

	host_mat[4]  = make_simpleComplexFloat( 0.058992f, -0.682839f);
	host_mat[5]  = make_simpleComplexFloat(-0.24367f,   0.851119f);
	host_mat[6]  = make_simpleComplexFloat( 1.63444f,   1.08207f);
	host_mat[7]  = make_simpleComplexFloat( 0.471581f, -0.283695f);

	host_mat[8]  = make_simpleComplexFloat( 0.427102f,  1.99722f);
	host_mat[9]  = make_simpleComplexFloat( 0.0900991f,-2.27228f);
	host_mat[10] = make_simpleComplexFloat( 0.294868f,  0.495371f);
	host_mat[11] = make_simpleComplexFloat(-1.07317f,  -0.438353f);

	host_mat[12] = make_simpleComplexFloat( 0.782054f, -1.08824f);
	host_mat[13] = make_simpleComplexFloat( 0.261584f, -0.122245f);
	host_mat[14] = make_simpleComplexFloat( 0.921061f,  0.306364f);
	host_mat[15] = make_simpleComplexFloat( 0.456282f,  2.01316f);

	//printf("\n");
	//print_square_matrix(host_mat, vector_dim);
	col_and_rows_swap_square_matrix(host_mat, vector_dim);
	//printf("transpose\n");
	//print_square_matrix(host_mat, vector_dim);

	/*
	column and rows swapped
	m = [ -0.736007-0.486131im  0.058992-0.682839im   0.427102+1.99722im   0.782054-1.08824im ;
		-0.50265+0.303371im  -0.24367+0.851119im  0.0900991-2.27228im   0.261584-0.122245im ;
		-1.13945-0.492444im   1.63444+1.08207im    0.294868+0.495371im  0.921061+0.306364im ;
		0.521185+1.68789im   0.471581-0.283695im   -1.07317-0.438353im  0.456282+2.01316im ]
	*/

	/*
	host_mat[0]  = make_simpleComplexFloat(-0.736007f, -0.486131f);
	host_mat[1]  = make_simpleComplexFloat( 0.058992f, -0.682839f);
	host_mat[2]  = make_simpleComplexFloat( 0.427102f,  1.99722f);
	host_mat[3]  = make_simpleComplexFloat( 0.782054f, -1.08824f);
	host_mat[4]  = make_simpleComplexFloat(-0.50265f,   0.303371f);
	host_mat[5]  = make_simpleComplexFloat(-0.24367f,   0.851119f);
	host_mat[6]  = make_simpleComplexFloat( 0.0900991f,-2.27228f);
	host_mat[7]  = make_simpleComplexFloat( 0.261584f, -0.122245f);
	host_mat[8]  = make_simpleComplexFloat(-1.13945f,  -0.492444f);
	host_mat[9]  = make_simpleComplexFloat( 1.63444f,   1.08207f);
	host_mat[10] = make_simpleComplexFloat( 0.294868f,  0.495371f);
	host_mat[11] = make_simpleComplexFloat( 0.921061f,  0.306364f);
	host_mat[12] = make_simpleComplexFloat( 0.521185f,  1.68789f);
	host_mat[13] = make_simpleComplexFloat( 0.471581f, -0.283695f);
	host_mat[14] = make_simpleComplexFloat(-1.07317f,  -0.438353f);
	host_mat[15] = make_simpleComplexFloat( 0.456282f,  2.01316f);
	*/

	//printf("vectors_set:\n");
	// int k;
	for (i = 0; i < sampling_points; i++)
	{
		/*
		for (k = 0; k < vector_dim; k++)
		{
			host_vectors_set[k + (i*vector_dim)] = 0.0f;
			host_tmpvs[k + (i*vector_dim)] = 0.0f;
		}
		*/
		host_pts[i].x = 1.0f;
		host_pts[i].y = 0.0f;
	}

 	int error_level = sampling_routine_simpleComplex_float(host_mat, host_pts, vector_dim, sampling_points, 0);

	printf("pts=[\n");
	for (i = 0; i < sampling_points; i++)
	{
		printf("%f+%fim ;\n", host_pts[i].x, host_pts[i].y);
	}
	printf("]\n");

	free((void*)host_pts);
	free((void*)host_mat);
}

int main(int argc, char *argv[])
{
	init_enviroment();

	sampling_routine_test_simpleComplex_float();

	/*
	simpleComplexFloat a, b, c;

	a = make_simpleComplexFloat(1.0f, 0.5f);
	b = make_simpleComplexFloat(0.25f, 0.75f);
	c = make_simpleComplexFloat(0.0f, 0.0f);

	c = a + b;
	*/

	finalize_environment();

	return 0;
}
