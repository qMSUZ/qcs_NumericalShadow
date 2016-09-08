#ifndef __samplingroutine_h__
#define __samplingroutine_h__

#include "cplxNum.h"


extern "C" void col_and_rows_swap_square_matrix(simpleComplexFloat *hm_in, const int lead_dim);
extern "C" void print_square_matrix(const simpleComplexFloat *hm, const int lead_dim);

extern "C" int sampling_routine_simpleComplex_float(simpleComplexFloat *host_mat,
													simpleComplexFloat *host_pts,
													const int vector_dim, const int sampling_points,
													const int sampling_type);


extern "C" int init_enviroment();
extern "C" int finalize_environment();

#endif // __samplingroutine_h__
