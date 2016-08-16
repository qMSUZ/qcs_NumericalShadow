#ifndef __samplingroutine_h__
#define __samplingroutine_h__

#include "cplxNum.cuh"


extern "C" void transpose_square_matrix(simpleComplex<float> *hm, const int lead_dim);

extern "C" int sampling_routine_simpleComplex_float(simpleComplex<float> *host_mat,
	simpleComplex<float> *host_pts,
	const int vector_dim, const int sampling_points, const int sampling_type);


extern "C" int init_enviroment();
extern "C" int finalize_environment();

#endif // __samplingroutine_h__
