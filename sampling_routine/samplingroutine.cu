#include <cstdio>
#include <ctime>


#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>


#include "cplxNum.cuh"
#include "cplxNum.cu"
#include "samplingroutine.cuh"

curandState *devStates;
//curandStateMRG32k3a  *devStates;

__global__ void setup_kernel_for_curand(curandState *state, unsigned long long seed)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	curand_init(seed, id, 0, &state[id]);
} /* end of kernel: setup_kernel_for_curand */

__global__ void kernel_random_quantum_state_simpleComplex_float(curandState *state, simpleComplex<float> *vector_state, float *rslt, const unsigned int N)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int stateID = blockIdx.x * blockDim.x + threadIdx.x;

	curandState localState = state[stateID];
	float sumTmp = 0.0f;
	float t = 0.0f;

	float rngval = (float)0.0;

	while(i < N)
	{
		rngval = (float)curand_uniform_double(&localState);
		vector_state[i].re = rngval;

		rngval = (float)curand_uniform_double(&localState);
		vector_state[i].im = rngval;

		t = simpleComplexMod( vector_state[i] );
		sumTmp += (t * t);

		i += blockDim.x*gridDim.x;
	}

	atomicAdd(rslt, sumTmp);
	
	state[stateID] = localState;
} /* end of kernel: kernel_random_quantum_state_simpleComplex_float */

__global__ void kernel_scale_vector_simpleComplex_float(simpleComplex<float> *a, float scale_val, const unsigned int N)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	while (i < N) {

		a[i] = a[i] * scale_val;

		i += gridDim.x * blockDim.x;
	}

	__syncthreads();
} /* end of kernel: kernel_scale_vector_simpleComplex_float */


/*
 * kernel: matvec_kernel_simpleComplex_float 

 Very important assumption:

	matrix is stored in column-major order !!!!

	to avoid problem with un-coalesced global memory access.

	More routines for matrix-vector multiplication can be found at:
	http://stackoverflow.com/questions/26417475/matrix-vector-multiplication-in-cuda-benchmarking-performance

*/

template<const unsigned int blk>
__global__ void matvec_kernel_simpleComplex_float(const simpleComplex<float> *dA,
							  const simpleComplex<float> *dx,
							  simpleComplex<float> *dy, const unsigned int nRows, const unsigned int nx)
{
	const unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
	const unsigned int hor_blocks = (nx + blk - 1) / blk;

	__shared__ simpleComplex<float> x_shared[blk];

	register simpleComplex<float> y_val;
	y_val.re = 0.0f;
	y_val.im = 0.0f;

#pragma unroll
	for (unsigned int m = 0; m < hor_blocks; ++m) {

		if ((m * blk + threadIdx.x) < nx)
		{
			x_shared[threadIdx.x] = dx[threadIdx.x + m * blk];

		}
		else
		{
			x_shared[threadIdx.x].re = 0.0f;
			x_shared[threadIdx.x].im = 0.0f;
		}

		__syncthreads();

#pragma unroll
		for (unsigned int e = 0; e < blk; ++e)
		{
			y_val = y_val + (dA[tid + (e + blk * m) * nRows] * x_shared[e]);
		}

		__syncthreads();
	}

	if (tid < nRows)
		dy[tid] = y_val;
} /* end of kernel: matvec_kernel_simpleComplex_float */


template<const size_t BLOCK_SIZE_DP>
__global__ void dot_product_simpleComplex_float(const simpleComplex<float>*  v1, const simpleComplex<float>*  v2, simpleComplex<float>*  rslt, int N)
{
	__shared__ simpleComplex<float> cache[BLOCK_SIZE_DP];
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	cache[threadIdx.x].re = 0.0f;
	cache[threadIdx.x].im = 0.0f;

	while (i < N) {
		cache[threadIdx.x] = cache[threadIdx.x] + (simpleComplexAdj(v1[i]) * v2[i]);
		i += gridDim.x * blockDim.x;
	}
	__syncthreads();

	i = BLOCK_SIZE_DP / 2;
	while (i > 0) {
		if (threadIdx.x < i) cache[threadIdx.x] = cache[threadIdx.x] + cache[threadIdx.x + i];
		__syncthreads();
		i /= 2;
	}

	if (threadIdx.x == 0)
	{
		atomicAdd(&rslt->re, cache[0].re);
		atomicAdd(&rslt->im, cache[0].im);
	}
}  /* end of kernel: dot_product_simpleComplex_float */


extern "C" void transpose_square_matrix(simpleComplex<float> *hm, const int lead_dim)
{
	int i, j;
	simpleComplex<float> a, b;

	for (i = 0; i < lead_dim; i++)
	{
		for (j = i; j < lead_dim; j++)
		{
			//printf("from (%d,%d)=%d to (%d,%d)=%d\n", i,j,i+(j*lead_dim), j,i,j + (i*lead_dim));
			a.re = hm[i+(j*lead_dim)].re;
			a.im = hm[i+(j*lead_dim)].im;

			//qcs_conj_complex(&a, &a);
			b.re = hm[j+(i*lead_dim)].re;
			b.im = hm[j+(i*lead_dim)].im;

			//qcs_conj_complex(&b, &b);
			//qcs_set_cell_at_matrix_complex(a_in, i, j, &b);
			hm[i + (j*lead_dim)].re = b.re;
			hm[i + (j*lead_dim)].im = b.im;

			//qcs_set_cell_at_matrix_complex(a_in, j, i, &a);
			hm[j + (i*lead_dim)].re = a.re;
			hm[j + (i*lead_dim)].im = a.im;
		}
	}
}

extern "C" int sampling_routine_simpleComplex_float(simpleComplex<float> *host_mat,
	simpleComplex<float> *host_pts,
	const int vector_dim, const int sampling_points, const int sampling_type)
{
	cudaError_t cudaStatus = cudaSuccess;
	int i, error_level = 0;

	unsigned long long seedtime = time(0);

	float rslt = 0.0f;
	simpleComplex<float> *dev_mat = NULL;
	simpleComplex<float> **dev_vectors_set;
	simpleComplex<float> **dev_tmpvs;
	float **dev_nrmval;
	simpleComplex<float> **dev_pts;

	dev_vectors_set = (simpleComplex<float> **)malloc(sizeof(simpleComplex<float>*) * sampling_points);
	dev_tmpvs = (simpleComplex<float> **)malloc(sizeof(simpleComplex<float>*) * sampling_points);
	dev_nrmval = (float **)malloc(sizeof(float*) * sampling_points);
	dev_pts = (simpleComplex<float> **)malloc(sizeof(simpleComplex<float>*) * sampling_points);

	curandState *devStates;

	cudaStatus = cudaMalloc((void **)&devStates, 64*16 * sizeof(curandState));
	if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

	cudaStatus = cudaMalloc((void **)&dev_mat, vector_dim * vector_dim * sizeof(simpleComplex<float>));
	if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

	for (i = 0; i < sampling_points; i++)
	{
		cudaStatus = cudaMalloc((void **)&dev_vectors_set[i], vector_dim *  sizeof(simpleComplex<float>));
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

		cudaStatus = cudaMalloc((void **)&dev_tmpvs[i], vector_dim * sizeof(simpleComplex<float>));
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

		cudaStatus = cudaMalloc((void **)&dev_nrmval[i], sizeof(float));
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

		cudaStatus = cudaMalloc((void **)&dev_pts[i],  sizeof(simpleComplex<float>));
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }
	}

	cudaStatus = cudaMemcpy(dev_mat, host_mat, vector_dim * vector_dim * sizeof(simpleComplex<float>), cudaMemcpyHostToDevice);
	cudaStatus = cudaDeviceSynchronize();

	setup_kernel_for_curand << <64, 16 >> >(devStates, seedtime);
	cudaStatus = cudaDeviceSynchronize();

	for (i = 0; i < sampling_points; i++)
	{
		kernel_random_quantum_state_simpleComplex_float << <64, 16, 16*sizeof(simpleComplex<float>) >> >(devStates, dev_vectors_set[i], dev_nrmval[i], vector_dim);
		cudaStatus = cudaDeviceSynchronize();

		cudaStatus = cudaMemcpy(&rslt, dev_nrmval[i], sizeof(float), cudaMemcpyDeviceToHost);
		cudaStatus = cudaDeviceSynchronize();

		//printf("nrmval[%d] %f\n", i, rslt);

		kernel_scale_vector_simpleComplex_float << <64, 16 >> >(dev_vectors_set[i], 1.0f/sqrt(rslt), vector_dim);
		cudaStatus = cudaDeviceSynchronize();

	}

	simpleComplex<float> *host_vectors_set = new simpleComplex<float>[vector_dim];
	float val=0.0f, nrmval = 0.0f;
	int k;
	for (i = 0; i < sampling_points; i++)
	{
		memset(host_vectors_set, 0, vector_dim * sizeof(simpleComplex<float>));
		cudaStatus = cudaMemcpy(host_vectors_set, dev_vectors_set[i], vector_dim * sizeof(simpleComplex<float>), cudaMemcpyDeviceToHost);
		cudaStatus = cudaDeviceSynchronize();

		/*
		printf("v%d=[", i); nrmval = 0.0f;
		for (k = 0; k < vector_dim; k++) {
			val = simpleComplexMod(host_vectors_set[k]);
			printf(" %f+%fim ;", host_vectors_set[k].re, host_vectors_set[k].im);
			nrmval += val*val;
		}
		printf("] nrm = %f\n", sqrtf(nrmval) );
		*/
	}
	delete[] host_vectors_set;

	const int blk = 16;

	dim3 dim_grid((vector_dim + blk - 1) / blk);
	dim3 dim_block(blk);

	for (i = 0; i < sampling_points; i++)
	{
		matvec_kernel_simpleComplex_float<blk> << <dim_grid, dim_block, sizeof(simpleComplex<float>)*blk >> >(dev_mat, dev_vectors_set[i], dev_tmpvs[i], vector_dim, vector_dim);
	}
	cudaDeviceSynchronize();

	for (i = 0; i < sampling_points; i++)
	{
		dot_product_simpleComplex_float<16> << <64, 16, sizeof(simpleComplex<float>)*16 >> >( dev_vectors_set[i], dev_tmpvs[i], dev_pts[i], vector_dim);
	}
	cudaDeviceSynchronize();

	for (i = 0; i < sampling_points; i++)
	{
		cudaStatus = cudaMemcpy(host_pts+i, dev_pts[i], sizeof(simpleComplex<float>), cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
	}

	cudaFree(dev_pts);
	
	for (i = 0; i < sampling_points; i++)
	{
		cudaFree(dev_nrmval[i]);
		cudaFree(dev_tmpvs[i]);
		cudaFree(dev_vectors_set[i]);
		cudaFree(dev_pts[i]);
	}
	cudaFree(dev_mat);
	cudaFree(devStates);

	free((void*)dev_pts);
	free((void*)dev_vectors_set);
	free((void*)dev_tmpvs);
	free((void*)dev_nrmval);

	return error_level;

} /* end of host function extern "C": sampling_routine_simpleComplex_float */


extern "C" int init_enviroment()
{
	printf("Sampling library is starting ...\n\n");
	printf("(+) CUDA ART version is %d.%d.\n", CUDART_VERSION / 1000, CUDART_VERSION % 1000);
	printf("(i) CUDA Device Query (Runtime API) version (CUDART static linking)\n");

	cudaSetDevice(0);

	return 0;

}
extern "C" int finalize_environment()
{
	printf("... sampling library is finalizing. \n\n");

	return 0;
}

