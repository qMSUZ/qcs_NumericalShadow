#ifndef __numrange_h__
#define __numrange_h__

#include "cublas_v2.h"
#include "magma_v2.h"
#include "magmablas_v2.h"
#include "magma_lapack.h"

#ifndef M_PI_F
#define M_PI_F       3.14159265358979323846f
#endif

#ifndef M_PI_D
#define M_PI_D       3.14159265358979323846
#endif

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

extern "C" int trival_method(magmaFloatComplex *M, magma_int_t N);
extern "C" int get_value();
extern "C" void set_value(int a);

extern "C" int init_enviroment();
extern "C" int finalize_environment();

extern "C" int calc_bounding_box(magmaFloatComplex *M, magma_int_t M_lead_dim, float *wReEig, float *wImEig);
extern "C" int calc_numerical_range(magmaFloatComplex *M, magma_int_t M_lead_dim, float _from, float _step, magma_int_t _steps, magmaFloatComplex *pts);

#endif // __numrange_h__
