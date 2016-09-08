#include <cstdio>
#include <cstring>
#include <ctime>


#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>


#include "cplxNum.h"
#include "samplingroutine.h"

curandState *devStates;
//curandStateMRG32k3a  *devStates;

__global__ void setup_kernel_for_curand(curandState *state, unsigned long long seed)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	//printf("setup_kernel_for_curand: id %d\n", id);

	curand_init(seed, id, 0, &state[id]);

} /* end of kernel: setup_kernel_for_curand */

__global__ void kernel_random_quantum_state_simpleComplex_float(curandState *state, simpleComplexFloat *vector_state, float *rslt, const unsigned int N)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	const int stateID = blockIdx.x * blockDim.x + threadIdx.x;

	curandState localState = state[stateID];

	simpleComplexFloat v;
	float sumTmp = 0.0f;
	float t = 0.0f;

	float rngval = (float)0.0;

	if (i == 0)
		*rslt = 0.0f;

	__syncthreads();

	while(i < N)
	{
		rngval = (float)curand_uniform_double(&localState);
		v.x = rngval;

		rngval = (float)curand_uniform_double(&localState);
		v.y = rngval;

		t = simpleComplexMod( v );
		//t = sqrtf(vector_state[i].x * vector_state[i].x + vector_state[i].y*vector_state[i].y);
		sumTmp += (t * t);

		vector_state[i].x = v.x;
		vector_state[i].y = v.y;

		//printf("i=%d, sumTmp=%f\n", i, sumTmp);
		//printf("i=%d, %f+%f*im ;\n", i, vector_state[i].x, vector_state[i].y);

		i += blockDim.x*gridDim.x;
	}
	state[stateID] = localState;

	atomicAdd(rslt, sumTmp);

} /* end of kernel: kernel_random_quantum_state_simpleComplex_float */

__global__ void kernel_scale_vector_simpleComplex_float( simpleComplexFloat *a, const float scale_val, const unsigned int N)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	while (i < N) {

		a[i] = a[i] * scale_val;

		i += blockDim.x*gridDim.x;
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
__global__ void matvec_kernel_simpleComplex_float(const  simpleComplexFloat *dA,
							  const  simpleComplexFloat *dx,
							   simpleComplexFloat *dy, const unsigned int nRows, const unsigned int nx)
{
	const unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
	const unsigned int hor_blocks = (nx + blk - 1) / blk;

	__shared__  simpleComplexFloat x_shared[blk];

	register  simpleComplexFloat y_val;
	y_val.x = 0.0f;
	y_val.y = 0.0f;

#pragma unroll
	for (unsigned int m = 0; m < hor_blocks; ++m) {

		if ((m * blk + threadIdx.x) < nx)
		{
			x_shared[threadIdx.x] = dx[threadIdx.x + m * blk];

		}
		else
		{
			x_shared[threadIdx.x].x = 0.0f;
			x_shared[threadIdx.x].y = 0.0f;
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
	{
		dy[tid] = y_val;
	}
} /* end of kernel: matvec_kernel_simpleComplex_float */


template<const size_t BLOCK_SIZE_DP>
__global__ void dot_product_simpleComplex_float(const  simpleComplexFloat *v1, const  simpleComplexFloat *v2,  simpleComplexFloat *rslt, int N)
{
	__shared__  simpleComplexFloat cache[BLOCK_SIZE_DP];
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	cache[threadIdx.x].x = 0.0f;
	cache[threadIdx.x].y = 0.0f;

	if (i == 0)
	{
		rslt->x = 0.0f;
		rslt->y = 0.0f;
	}

	__syncthreads();
	

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
		atomicAdd(&rslt->x, cache[0].x);
		atomicAdd(&rslt->y, cache[0].y);
	}
}  /* end of kernel: dot_product_simpleComplex_float */


/* transpose without adjoint */
extern "C" void col_and_rows_swap_square_matrix( simpleComplexFloat *hm_in, const int lead_dim)
{
	int i, j;
	 simpleComplexFloat a, b;

#define _addrElem(x,y) (x*lead_dim)+y

	for (i = 0; i < lead_dim; i++) // row
	{
		for (j = i; j < lead_dim; j++) // col
		{
			//printf("from (%d,%d)=%d to (%d,%d)=%d\n", i,j,i+(j*lead_dim), j,i,j + (i*lead_dim));
			a.x = hm_in[_addrElem(i,j)].x;
			a.y = hm_in[_addrElem(i,j)].y;

			b.x = hm_in[_addrElem(j,i)].x;
			b.y = hm_in[_addrElem(j,i)].y;

			hm_in[_addrElem(i,j)].x = b.x;
			hm_in[_addrElem(i,j)].y = b.y;

			hm_in[_addrElem(j,i)].x = a.x;
			hm_in[_addrElem(j,i)].y = a.y;


			//qcs_conj_complex(&a, &a);
			//b.x = hm[j+(i*lead_dim)].x;
			//b.y = hm[j+(i*lead_dim)].y;

			//qcs_conj_complex(&b, &b);
			//qcs_set_cell_at_matrix_complex(a_in, i, j, &b);
			//hm[i + (j*lead_dim)].x = b.x;
			//hm[i + (j*lead_dim)].x = b.y;

			//qcs_set_cell_at_matrix_complex(a_in, j, i, &a);
			//hm[j + (i*lead_dim)].x = a.x;
			//hm[j + (i*lead_dim)].y = a.y;
		}
	}
#undef _addrElem
}

extern "C" void print_square_matrix( const simpleComplexFloat *hm, const int lead_dim)
{
	int i, j;
	simpleComplexFloat a, b;

	printf("m = [");
	for (i = 0; i < lead_dim; i++) // row
	{
		for (j = 0; j < lead_dim; j++) // col
		{
			printf("%f+%fim ", hm[j+(i*lead_dim)].x, hm[j+(i*lead_dim)].y);
		}
		printf(";\n");
	}
	printf("]\n");
}


extern "C" int sampling_routine_simpleComplex_float( simpleComplexFloat *host_mat,
	 simpleComplexFloat *host_pts,
	const int vector_dim, const int sampling_points, const int sampling_type)
{
	cudaError_t cudaStatus = cudaSuccess;
	int i=0, k=0, error_level = 0;

	unsigned long long seedtime = time(0);

	simpleComplexFloat tmp_vector_set[4];

	float rslt = 0.0f;
	simpleComplexFloat *dev_mat = NULL;
	simpleComplexFloat *dev_vectors_set = NULL;
	simpleComplexFloat *dev_tmpvs = NULL;
	float *dev_nrmval = NULL;
	simpleComplexFloat *dev_pts = NULL;

	float val = 0.0f, nrmval = 0.0f;

	//dev_vectors_set = (simpleComplexFloat **)malloc(sizeof(simpleComplexFloat*) * sampling_points);
	//dev_tmpvs = (simpleComplexFloat **)malloc(sizeof(simpleComplexFloat*) * sampling_points);
	//dev_nrmval = (float **)malloc(sizeof(float*) * sampling_points);
	//dev_pts = ( simpleComplexFloat **)malloc(sizeof( simpleComplexFloat*) * sampling_points);

	curandState *devStates;

	cudaStatus = cudaMalloc(&devStates, 64*16*sizeof(curandState));
	if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

	cudaStatus = cudaMalloc((void **)&dev_mat, vector_dim * vector_dim * sizeof(simpleComplexFloat));
	if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

	cudaStatus = cudaMalloc((void **)&dev_vectors_set, vector_dim *  sizeof(simpleComplexFloat));
	if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

	cudaStatus = cudaMalloc((void **)&dev_tmpvs, vector_dim * sizeof(simpleComplexFloat));
	if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

	cudaStatus = cudaMalloc((void **)&dev_nrmval, sizeof(float));
	if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

	cudaStatus = cudaMalloc((void **)&dev_pts,  sizeof(simpleComplexFloat));
	if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaMalloc failed!"); return cudaStatus; }

	cudaStatus = cudaMemcpy(dev_mat, host_mat, vector_dim * vector_dim * sizeof(simpleComplexFloat), cudaMemcpyHostToDevice);
	cudaStatus = cudaDeviceSynchronize();

	setup_kernel_for_curand << <64, 16 >> >(devStates, seedtime);
	cudaStatus = cudaDeviceSynchronize();

	for (i = 0; i < sampling_points; i++)
	{
		kernel_random_quantum_state_simpleComplex_float << <64, 16 >> >(devStates, dev_vectors_set, dev_nrmval, vector_dim);
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "kernel_random_quantum_state_simpleComplex_float failed!"); return cudaStatus; }


		/// debug code begin
#if 0		
		cudaStatus = cudaMemcpy(&tmp_vector_set, dev_vectors_set, vector_dim * sizeof(simpleComplexFloat), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "tmp_vector_set: cudaMemcpy failed!"); return cudaStatus; }
		cudaStatus = cudaDeviceSynchronize();
		printf("tv%d=[", i);
		for (k = 0; k < vector_dim; k++)
		{
			printf("%f+%fim; ", tmp_vector_set[k].x, tmp_vector_set[k].y);
		}
		printf("]; \n");
#endif
		/// debug code end

		cudaStatus = cudaMemcpy(&rslt, dev_nrmval, sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "rslt: cudaMemcpy failed!"); return cudaStatus; }
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaDeviceSynchronize failed!"); return cudaStatus; }

		/// debug code begin
#if 0		
		printf("nrmval%d=%f;\n", i, sqrtf(rslt));
#endif
		/// debug code end

		kernel_scale_vector_simpleComplex_float << <64, 16 >> >(dev_vectors_set, 1.0f / sqrt(rslt), vector_dim);
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "kernel_scale_vector_simpleComplex_float failed!"); return cudaStatus; }
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) { fprintf(stderr, "cudaDeviceSynchronize failed!"); return cudaStatus; }


		// simpleComplexFloat *host_vectors_set = new  simpleComplexFloat[vector_dim];
		//for (i = 0; i < sampling_points; i++)
		
		/// debug code begin
		memset(&tmp_vector_set, 0, vector_dim * sizeof(simpleComplexFloat));
		cudaStatus = cudaMemcpy(&tmp_vector_set, dev_vectors_set, vector_dim * sizeof(simpleComplexFloat), cudaMemcpyDeviceToHost);
		cudaStatus = cudaDeviceSynchronize();

		printf("v%d=[", i); 
		nrmval = 0.0f;
		for (k = 0; k < vector_dim; k++) {
			val = simpleComplexMod(tmp_vector_set[k]);
			//val = sqrtf(host_vectors_set[k].x * host_vectors_set[k].x + host_vectors_set[k].y*host_vectors_set[k].y);
			printf(" %f+%fim ;", tmp_vector_set[k].x, tmp_vector_set[k].y);
			nrmval += val*val;
		}
		printf("]; nrmv%d = %f;\n", i, sqrtf(nrmval));
		/// debug code end

		//delete[] host_vectors_set;

		const int blk = 16;
		dim3 dim_grid((vector_dim + blk - 1) / blk);
		dim3 dim_block(blk);

		matvec_kernel_simpleComplex_float<blk> << <dim_grid, dim_block, sizeof(simpleComplexFloat)*blk >> >(dev_mat, dev_vectors_set, dev_tmpvs, vector_dim, vector_dim);
		cudaDeviceSynchronize();

		dot_product_simpleComplex_float<16> << <64, 16, sizeof(simpleComplexFloat) * 16 >> >(dev_vectors_set, dev_tmpvs, dev_pts, vector_dim);
		cudaDeviceSynchronize();

		cudaStatus = cudaMemcpy(&host_pts[i], dev_pts, sizeof(simpleComplexFloat), cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
	} // for (i = 0; i < sampling_points; i++)

	
	cudaFree(dev_nrmval);
	cudaFree(dev_tmpvs);
	cudaFree(dev_vectors_set);
	cudaFree(dev_pts);

	cudaFree(dev_mat);
	cudaFree(devStates);

	//free((void*)dev_pts);
	//free((void*)dev_vectors_set);
	//free((void*)dev_tmpvs);
	//free((void*)dev_nrmval);

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
