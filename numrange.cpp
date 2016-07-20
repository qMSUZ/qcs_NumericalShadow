#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <complex>


#include <algorithm>

#include "cublas_v2.h"
#include "magma_v2.h"
#include "magmablas_v2.h"
#include "magma_lapack.h"

#include "numrange.h"

static magma_int_t dev = 0;
static magma_queue_t queue = NULL;

static magma_int_t value = 0;

extern "C" int get_value()
{
	return value;
}

extern "C" void set_value(int a)
{
	 value = a;
}

extern "C" int trival_method(magmaFloatComplex *M, magma_int_t N)
{
	int i;

	for(i=0;i<N;i++)
	{
		printf("%4d: %f %f\n", i, M[i].x, M[i].y);
	}

	return 0;
}

extern "C" int init_enviroment()
{
	magma_init();
	magma_print_environment();

	magma_queue_create(dev, &queue);


	return 0;
}


extern "C" int finalize_environment()
{
	magma_queue_destroy(queue);

	magma_finalize();

	return 0;
}

extern "C" int calc_bounding_box(magmaFloatComplex *M, magma_int_t M_lead_dim, float *wReEig, float *wImEig)
{
	magma_int_t rslt = 0;

	//magmaFloatComplex *AT = nullptr;
	magmaFloatComplex *dA = nullptr, *dAT = nullptr,
		*dreA = nullptr, *dimA = nullptr;

	float *dreEig = nullptr;
	float *dimEig = nullptr;

	//magma_int_t *ipiv = NULL;
	magma_int_t lda = M_lead_dim;
	//magma_int_t ldx = lda;
	magma_int_t info = 0;

	magma_int_t nb = 0;

	//magma_vec_t jobvl;
	//magma_vec_t jobvr;

	magmaFloatComplex *work = nullptr;
	magma_int_t  lwork = 0;

	float *rwork = nullptr;
	magma_int_t lrwork = 0;

	magma_int_t *iwork = nullptr;
	magma_int_t liwork = 0;


	nb = magma_get_cgehrd_nb( M_lead_dim );


	lwork = 2 * (M_lead_dim + M_lead_dim*nb); // MagmaNoVec
	//lwork = 2 * max(M_lead_dim + M_lead_dim*nb, 2*M_lead_dim + M_lead_dim*M_lead_dim); // MagmaVec

	lrwork = M_lead_dim; // MagmaNoVec
	//lrwork = 1 + 5 * M_lead_dim + 2*M_lead_dim*M_lead_dim; // MagmaVec

	liwork = 1; // MagmaNoVec
	//liwork = 3 + 5*M_lead_dim; // MagmaVec

	magma_imalloc_cpu(&iwork, liwork);

	magma_smalloc_cpu(&rwork, lrwork);

	//magma_cmalloc_cpu(&A, lda*M_lead_dim);
	//magma_cmalloc_cpu(&AT, lda*M_lead_dim);

	//magma_smalloc_cpu(&reEig, M_lead_dim);
	//magma_smalloc_cpu(&imEig, M_lead_dim);


	magma_cmalloc_pinned(&dA, lda*M_lead_dim);
	magma_cmalloc_pinned(&dAT, lda*M_lead_dim);

	magma_cmalloc_pinned(&dreA, lda*M_lead_dim);
	magma_cmalloc_pinned(&dimA, lda*M_lead_dim);

	//magma_cmalloc_pinned(&VL, lda*M_lead_dim);
	//magma_cmalloc_pinned(&VR, lda*M_lead_dim);

	magma_cmalloc_pinned(&work, lwork);

	magma_smalloc_pinned(&dreEig, M_lead_dim);
	magma_smalloc_pinned(&dimEig, M_lead_dim);


	//matrix_fillzero(AT, M_lead_dim);

	//vector_fillzero(reEig, M_lead_dim);
	//vector_fillzero(imEig, M_lead_dim);

	//prepare_matrix_2(M);

	magma_csetmatrix(M_lead_dim, M_lead_dim, M, lda, dA, M_lead_dim, queue);
	//magma_csetmatrix(M_lead_dim, M_lead_dim, AT, lda, dAT, M_lead_dim, queue);

	//magma_ssetvector(M_lead_dim, wReEig, 1, dreEig, 1, queue);
	//magma_ssetvector(M_lead_dim, wImEig, 1, dimEig, 1, queue);

	//magma_cprint_gpu(M_lead_dim, M_lead_dim, dA, lda);

	// reA = ( (A + A')/2.0 )
	// A'
	magmablas_ctranspose(M_lead_dim, M_lead_dim, dA, M_lead_dim, dAT, M_lead_dim, queue);
	//magma_cprint_gpu(M_lead_dim, M_lead_dim, dAT, lda);

	// AT = A + A'
	magmablas_cgeadd(M_lead_dim, M_lead_dim, MAGMA_C_MAKE(1.0f, 0.0f), dA, M_lead_dim, dAT, M_lead_dim, queue);
	//magma_cprint_gpu(M_lead_dim, M_lead_dim, dAT, lda);
	// AT=AT*0.5
	magma_cscal(lda*M_lead_dim, MAGMA_C_MAKE(0.5f, 0.0f), dAT, 1, queue);
	//magma_cprint_gpu(M_lead_dim, M_lead_dim, dAT, lda);
	// reA = AT
	magma_ccopy(lda*M_lead_dim, dAT, 1, dreA, 1, queue);
	//magma_cprint_gpu(M_lead_dim, M_lead_dim, dreA, lda);
	magma_sync_wtime(queue);


	// imA = ( -1im*(A - A')/2.0 )
	// A'
	magmablas_ctranspose(M_lead_dim, M_lead_dim, dA, M_lead_dim, dAT, M_lead_dim, queue);
	//magma_cprint_gpu(M_lead_dim, M_lead_dim, dAT, lda);
	// AT = A + A'
	magmablas_cgeadd(M_lead_dim, M_lead_dim, MAGMA_C_MAKE(-1.0f, 0.0f), dAT, M_lead_dim, dA, M_lead_dim, queue);
	// A=A*-1j*0.5
	magma_cscal(lda*M_lead_dim, MAGMA_C_MAKE(0.0f, -0.5f), dA, 1, queue);
	// imA = A
	magma_ccopy(lda*M_lead_dim, dA, 1, dimA, 1, queue);
	magma_sync_wtime(queue);

	//magma_cprint_gpu(M_lead_dim, M_lead_dim, dreA, lda);
	//magma_cprint_gpu(M_lead_dim, M_lead_dim, dimA, lda);


	// reEig::Vector=eigvals(reA)
	rslt = magma_cheevd(MagmaNoVec, MagmaLower,
		M_lead_dim,
		dreA, lda,
		dreEig,
		work, lwork,
		rwork, lrwork,
		iwork, liwork,
		&info);

	// imEig::Vector=eigvals(imA)
	rslt = magma_cheevd(MagmaNoVec, MagmaLower,
		M_lead_dim,
		dimA, lda,
		dimEig,
		work, lwork,
		rwork, lrwork,
		iwork, liwork,
		&info);


	//magma_sprint_gpu(M_lead_dim, 1, dreEig, M_lead_dim);
	//magma_sprint_gpu(M_lead_dim, 1, dimEig, M_lead_dim);


	magma_sgetvector(M_lead_dim, dreEig, 1, wReEig, 1, queue);
	//magma_sync_wtime(queue);

	magma_sgetvector(M_lead_dim, dimEig, 1, wImEig, 1, queue);
	//magma_sync_wtime(queue);

	/*
	maxReIdx = magma_isamax(M_lead_dim, dreEig, 1, queue) - 1;
	minReIdx = magma_isamin(M_lead_dim, dreEig, 1, queue) - 1;

	maxImIdx = magma_isamax(M_lead_dim, dimEig, 1, queue) - 1;
	minImIdx = magma_isamin(M_lead_dim, dimEig, 1, queue) - 1;


	printf("max re idx = %d\nmin re idx = %d\n", maxReIdx, minReIdx);
	printf("%f %f\n", wReEig[maxReIdx], wReEig[minReIdx]);

	printf("max im idx = %d\nmin im idx = %d\n", maxImIdx, minImIdx);
	printf("%f %f\n", wImEig[maxImIdx], wImEig[minImIdx]);
	*/

	//printf("test wReEig: %f %f\n", wReEig[0], wReEig[1]);
	//printf("test wImEig: %f %f\n", wImEig[0], wImEig[1]);


	magma_free_cpu(iwork);
	magma_free_cpu(rwork);
	//magma_free_cpu(AT);

	magma_free_pinned(dA);
	magma_free_pinned(dAT);

	magma_free_pinned(dreA);
	magma_free_pinned(dimA);

	magma_free_pinned(work);

	magma_free_pinned(dreEig);
	magma_free_pinned(dimEig);

	return rslt;
}


extern "C" int calc_numerical_range(magmaFloatComplex *M, magma_int_t M_lead_dim, float _from, float _step, magma_int_t _steps, magmaFloatComplex *pts)
{
	magma_int_t idx = 0, rslt = 0;

	magmaFloatComplex p, scalar;
	std::complex<float> vtmp;

	float j;

	magmaFloatComplex *dA = nullptr;
	magmaFloatComplex *dAth = NULL, *dAthT = NULL,
				*dX = NULL, *dY = NULL;

	float *dE = NULL;
	//float *hE = NULL;


	//magma_int_t *ipiv = NULL;
	magma_int_t lda = M_lead_dim;
	//magma_int_t ldx = lda;
	magma_int_t info = 0;

	magma_int_t nb = 0;

	//magma_vec_t jobvl;
	//magma_vec_t jobvr;

	magmaFloatComplex *work = nullptr;
	magma_int_t  lwork = 0;

	float *rwork = nullptr;
	magma_int_t lrwork = 0;

	magma_int_t *iwork = nullptr;
	magma_int_t liwork = 0;

	nb = magma_get_cgehrd_nb( M_lead_dim );

	lwork = 2 * max(M_lead_dim + M_lead_dim*nb, 2 * M_lead_dim + M_lead_dim*M_lead_dim); // MagmaVec

	lrwork = 1 + 5 * M_lead_dim + 2 * M_lead_dim*M_lead_dim; // MagmaVec

	liwork = (3 + 5 * M_lead_dim); // MagmaVec

	magma_imalloc_cpu(&iwork, liwork);
	magma_smalloc_cpu(&rwork, lrwork);

	magma_cmalloc_pinned(&work, lwork);

	magma_cmalloc_pinned(&dA, lda*M_lead_dim);
	magma_cmalloc_pinned(&dAth, lda*M_lead_dim);
	magma_cmalloc_pinned(&dAthT, lda*M_lead_dim);

	magma_smalloc_pinned(&dE, M_lead_dim);
	//magma_smalloc_cpu(&hE, M_lead_dim);

	magma_cmalloc_pinned(&dX, M_lead_dim);
	magma_cmalloc_pinned(&dY, M_lead_dim);

	magma_csetmatrix(M_lead_dim, M_lead_dim, M, lda, dA, M_lead_dim, queue);

	// th=[0:resolution:2*pi]
	j = _from;
	for (idx = 0; idx < _steps; idx++)
	{
		//scalar = exp( 1im * -j);
		vtmp.real( 0.0f );
		vtmp.imag(  -j  );
		//vtmp = _FCbuild(0.0f, -j);
		//printf("vtmp = %f + i%f\n", vtmp._Val[0], vtmp._Val[1]);

		vtmp = exp(vtmp);
		scalar.x = vtmp.real();
		scalar.y = vtmp.imag();

		//printf("scalar = %f + i%f\n", scalar.x, scalar.y);

		magma_ccopy(lda * M_lead_dim, dA, 1, dAth, 1, queue);
		// Ath = exp(1im * -j) * As
		magma_cscal(lda * M_lead_dim, scalar, dAth, 1, queue);

		//magma_cprint_gpu(N, N, dA, lda);
		//magma_cprint_gpu(N, N, dAth, lda);

		// AthT = (Ath + Ath')
		magmablas_ctranspose_conj(M_lead_dim, M_lead_dim, dAth, M_lead_dim, dAthT, M_lead_dim, queue);
		magmablas_cgeadd(M_lead_dim, M_lead_dim, MAGMA_C_MAKE(1.0f, 0.0f), dAth, M_lead_dim, dAthT, M_lead_dim, queue);
		// AthT = AthT / 2
		magma_cscal(lda*M_lead_dim, MAGMA_C_MAKE(0.5f, 0.0f), dAthT, 1, queue);
		magma_sync_wtime(queue);

		//magma_cprint_gpu(M_lead_dim, M_lead_dim, dAthT, lda);

		// e, r = eig(AthT)
		rslt = magma_cheevd(MagmaVec, MagmaLower,
			M_lead_dim,
			dAthT, lda,
			dE,
			work, lwork,
			rwork, lrwork,
			iwork, liwork,
			&info);
		magma_sync_wtime(queue);

		//printf("magma_cheevd info=%d\n", info);

		//magma_cprint_gpu(M_lead_dim, M_lead_dim, dAthT, lda);
		//magma_sprint_gpu(M_lead_dim, 1, dE, M_lead_dim);

		//magma_sgetvector(M_lead_dim, dE, 1, hE, 1, queue);

		//printf("%f %f\n", hE[0], hE[1]);

		// p = r[:,s]' * A * r[:,s]
		// r = r[:,s]
		magma_ccopy(
			M_lead_dim,
			dAthT + (M_lead_dim*(M_lead_dim-1)), 1, // dAthT + (N), where (N) is a column offset
			dX, 1,
			queue);
		magma_sync_wtime(queue);

		//magma_cprint_gpu(M_lead_dim, 1, dX, M_lead_dim);

		// pp = A * r[:,s]
		magma_cgemv(MagmaNoTrans,
			M_lead_dim, M_lead_dim,
			MAGMA_C_MAKE(1.0f, 0.0f),
			dA, lda,
			dX, 1,
			MAGMA_C_MAKE(0.0f, 0.0f),
			dY, 1, queue);
		magma_sync_wtime(queue);

		//magma_cprint_gpu(M_lead_dim, 1, dY, M_lead_dim);

		// p = r' * pp
		p = magma_cdotc(M_lead_dim, dX, 1, dY, 1, queue);
		magma_sync_wtime(queue);

		pts[idx] = p;

		//printf("p = %f %fi\n", p.x, p.y);

		j += _step;
	} // end of for (idx = 0; idx < _steps; idx++)

	magma_free_pinned(dY);
	magma_free_pinned(dX);

	//magma_free_cpu(hE);
	magma_free_pinned(dE);

	magma_free_pinned(dAthT);
	magma_free_pinned(dAth);
	magma_free_pinned(dA);

	magma_free_pinned(work);

	magma_free_cpu(rwork);
	magma_free_cpu(iwork);
	//magma_free_cpu(w);
	//magma_free_cpu(A);

	return rslt;
}
