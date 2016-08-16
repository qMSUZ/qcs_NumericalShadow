#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ccomplex>


#include <algorithm>

#include "cublas_v2.h"
#include "magma_v2.h"
#include "magmablas_v2.h"
#include "magma_lapack.h"

#include "numrange.h"


void prepare_matrix_1(magmaFloatComplex *A, int lda)
{
#define A(i_, j_) A[ (i_) + (j_)*lda ]
	A(0, 0) = MAGMA_C_MAKE(0.806581f, 0.0f);
	A(0, 1) = MAGMA_C_MAKE(0.814736f, 0.0f);
	A(1, 0) = MAGMA_C_MAKE(0.679601f, 0.0f);
	A(1, 1) = MAGMA_C_MAKE(0.0725491f, 0.0f);
#undef A
}

/*
A = [ 0.5 0.25 ; 0.25 0.5 ];
*/


void prepare_matrix_2(magmaFloatComplex *A, int lda)
{
#define A(i_, j_) A[ (i_) + (j_)*lda ]
	A(0, 0) = MAGMA_C_MAKE(0.5f, 0.0f);
	A(0, 1) = MAGMA_C_MAKE(0.25f, 0.0f);
	A(1, 0) = MAGMA_C_MAKE(0.25f, 0.0f);
	A(1, 1) = MAGMA_C_MAKE(0.5f, 0.0f);
#undef A
}


void prepare_matrix_2b(magmaFloatComplex *A, int lda)
{
#define A(i_, j_) A[ (i_) + (j_)*lda ]
	A(0, 0) = MAGMA_C_MAKE( 0.231029f, 0.0f);
	A(0, 1) = MAGMA_C_MAKE(-1.540030f, 0.0f);
	A(1, 0) = MAGMA_C_MAKE( 0.864749f, 0.0f);
	A(1, 1) = MAGMA_C_MAKE( 0.857362f, 0.0f);
#undef A
}

/*
A=[-1.553230 1.433990 ;
0.231705 1.061290]
*/
void prepare_matrix_2c(magmaFloatComplex *A, int lda)
{
#define A(i_, j_) A[ (i_) + (j_)*lda ]
	A(0, 0) = MAGMA_C_MAKE(-1.553230f, 0.0f);
	A(0, 1) = MAGMA_C_MAKE( 1.433990f, 0.0f);
	A(1, 0) = MAGMA_C_MAKE( 0.231705f, 0.0f);
	A(1, 1) = MAGMA_C_MAKE( 1.061290f, 0.0f);
#undef A
}


void prepare_matrix_3(magmaFloatComplex *A, int lda)
{
#define A(i_, j_) A[ (i_) + (j_)*lda ]
	A(0, 0) = MAGMA_C_MAKE(-0.398287f, 0.0f);
	A(0, 1) = MAGMA_C_MAKE(-1.612910f, 0.0f);
	A(0, 2) = MAGMA_C_MAKE( 0.246223f, 0.0f);
	A(1, 0) = MAGMA_C_MAKE( 0.124418f, 0.0f);
	A(1, 1) = MAGMA_C_MAKE(-1.209970f, 0.0f);
	A(1, 2) = MAGMA_C_MAKE(-1.268770f, 0.0f);
	A(2, 0) = MAGMA_C_MAKE(-0.368598f, 0.0f);
	A(2, 1) = MAGMA_C_MAKE( 1.325440f, 0.0f);
	A(2, 2) = MAGMA_C_MAKE(-0.986160f, 0.0f);
#undef A
}

/*
A = [
2.659990  1.941440  0.856206 -2.092150;
1.848220 -0.943889  0.494220  0.549636;
0.735044 -0.524624 -0.819604  0.480188;
0.291891 -0.779152 -0.482638 -1.238260
]
*/

void prepare_matrix_4(magmaFloatComplex *A, int lda)
{
#define A(i_, j_) A[ (i_) + (j_)*lda ]
	A(0, 0) = MAGMA_C_MAKE( 2.659990f, 0.0f);
	A(0, 1) = MAGMA_C_MAKE( 1.941440f, 0.0f);
	A(0, 2) = MAGMA_C_MAKE( 0.856206f, 0.0f);
	A(0, 3) = MAGMA_C_MAKE(-2.092150f, 0.0f);
	A(1, 0) = MAGMA_C_MAKE( 1.848220f, 0.0f);
	A(1, 1) = MAGMA_C_MAKE(-0.943889f, 0.0f);
	A(1, 2) = MAGMA_C_MAKE( 0.494220f, 0.0f);
	A(1, 3) = MAGMA_C_MAKE( 0.549636f, 0.0f);
	A(2, 0) = MAGMA_C_MAKE( 0.735044f, 0.0f);
	A(2, 1) = MAGMA_C_MAKE(-0.524624f, 0.0f);
	A(2, 2) = MAGMA_C_MAKE(-0.819604f, 0.0f);
	A(2, 3) = MAGMA_C_MAKE( 0.480188f, 0.0f);
	A(3, 0) = MAGMA_C_MAKE( 0.291891f, 0.0f);
	A(3, 1) = MAGMA_C_MAKE(-0.779152f, 0.0f);
	A(3, 2) = MAGMA_C_MAKE(-0.482638f, 0.0f);
	A(3, 3) = MAGMA_C_MAKE(-1.238260f, 0.0f);
#undef A
}


void matrix_fillzero(magmaFloatComplex *A, int lda)
{
#define A(i_, j_) A[ (i_) + (j_)*lda ]
	A(0, 0) = MAGMA_C_MAKE(0.0f, 0.0f);
	A(0, 1) = MAGMA_C_MAKE(0.0f, 0.0f);
	A(1, 0) = MAGMA_C_MAKE(0.0f, 0.0f);
	A(1, 1) = MAGMA_C_MAKE(0.0f, 0.0f);
#undef A
}

void matrix_fillzero(float *A, int lda)
{
#define A(i_, j_) A[ (i_) + (j_)*lda ]
	A(0, 0) = 0.0f;
	A(0, 1) = 0.0f;
	A(1, 0) = 0.0f;
	A(1, 1) = 0.0f;
#undef A
}


void vector_fillzero(magmaFloatComplex *A, magma_int_t N)
{
	int i;

	for (i = 0; i < N; i++)
	{
		A[i] = MAGMA_C_MAKE(0.0f, 0.0f);
	}
}

void vector_fillzero(float *A, magma_int_t N)
{
	int i;

	for (i = 0; i < N; i++)
	{
		A[i] = 0.0f;
	}
}

void calc_bounding_box_test()
{
	int rslt, N = 4;

	magmaFloatComplex *A = nullptr;
	float *reEig = nullptr;
	float *imEig = nullptr;

	magma_cmalloc_cpu(&A, N*N);
	magma_smalloc_cpu(&reEig, N);
	magma_smalloc_cpu(&imEig, N);

	vector_fillzero(reEig, N);
	vector_fillzero(imEig, N);

	matrix_fillzero(A, N);

	//prepare_matrix_2(A, N);
	// result for 2:
	//   0.75
	//   0.25
	//   0.00
	//   0.00

	//prepare_matrix_2b(A, N);
	// result for 2b:
	//	 0.0836791
	//	 1.00471
	//	-1.20239
	//	 1.20239

	//prepare_matrix_3(A, N);
	// result for 3:
	//	-1.651910
	//	 0.047931
	//	-1.591090
	//	 1.591090

	prepare_matrix_4(A, N);
	// result for 3:
	//	-1.651910
	//	 0.047931
	//	-1.591090
	//	 1.591090


	rslt = calc_bounding_box(A, N, reEig, imEig);

	printf("rslt = %d\n", rslt);

	printf("re values: %f %f %f %f\n", reEig[0], reEig[1], reEig[2], reEig[3]);
	printf("im values: %f %f %f %f\n", imEig[0], imEig[1], imEig[2], imEig[3]);


	magma_free_cpu(imEig);
	magma_free_cpu(reEig);
	magma_free_cpu(A);
}

void calc_numerical_range_test()
{
	const int N = 2;

	magmaFloatComplex *A = nullptr;
	magmaFloatComplex *pts = nullptr;

	magma_cmalloc_cpu(&A,  N*N );

	matrix_fillzero(A, N);


	float from_val, to_val, step_val;
	int i, rslt, points;

	from_val = 0.0f;
	to_val = 2.0f * M_PI_F;
	step_val = 0.01f;

	points = 1 + (int)((to_val - from_val) / step_val);

	magma_cmalloc_cpu(&pts,points);

	vector_fillzero(pts, points);

	printf("number of points = %d\n", points);

	prepare_matrix_2(A, N);
	//prepare_matrix_2c(A, N);
	//prepare_matrix_4(A, N);

	rslt = calc_numerical_range(A, N, from_val, step_val, points, pts);

	for (i = 0; i < points; i++)
	{
		printf("%4d: %f+i%f ", i, pts[i].x, pts[i].y);
		if (((i+1) % 4) == 0) printf("\n");
	}

	magma_free_cpu(pts);
	magma_free_cpu(A);
}

int main(int argc, char *argv[])
{	
	int r = 0;

	r = init_enviroment();

	//calc_bounding_box_test();
	//calc_bounding_box_gpu_test_old();

	//calc_numerical_range();
	calc_numerical_range_test();

	r = finalize_environment();

	return r;
}
