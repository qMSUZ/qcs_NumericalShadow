#include "cplxNum.h"

#if 0
template <typename T>
__host__ __device__ struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a.re + b.re;
	t.im = a.im + b.im;

	return t;
}

template <typename T>
__host__ __device__ struct simpleComplex<T> operator-(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a.re - b.re;
	t.im = a.im - b.im;

	return t;
}

template <typename T>
__host__ __device__ struct simpleComplex<T> operator+(const T &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a + b.re;
	t.im = b.im;

	return t;
}

template <typename T>
__host__ __device__  struct simpleComplex<T> operator-(const T &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a - b.re;
	t.im = b.im;

	return t;
}

template <typename T>
__host__ __device__  struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const T &b)
{
	struct simpleComplex<T> t;

	t.re = a.re + b;
	t.im = a.im ;

	return t;
}

template <typename T>
__host__ __device__  struct simpleComplex<T> operator-(const struct simpleComplex<T> &a, const T &b)
{
	struct simpleComplex<T> t;

	t.re = a.re - b;
	t.im = a.im ;

	return t;
}

template <typename T>
__host__ __device__  struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = (a.re * b.re) - (a.im * b.im);
	t.im = (a.re * b.im) + (a.im * b.re);

	return t;
}

template <typename T>
__host__ __device__  struct simpleComplex<T> operator*(const T &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = (a * b.re);
	t.im = (a * b.im);

	return t;
}

template <typename T>
__host__ __device__  struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const T &b)
{
	struct simpleComplex<T> t;

	t.re = (a.re * b);
	t.im = (a.im * b);

	return t;
}


template <typename T>
__host__ __device__  struct simpleComplex<T> operator/(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	T s =  (b.re * b.re) + (b.im * b.im);

	t.re = ((a.re * b.re) + (a.im * b.im)) / s;
	t.im = ((a.im * b.re) - (a.re * b.im)) / s;

	return t;
}

template <typename T>
__host__ __device__  struct simpleComplex<T> operator/(const T &_a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t, a;

	a.re = a;
	a.im = (T)0.0;

	T s =  (b.re * b.re) + (b.im * b.im);

	t.re = ((a.re * b.re) + (a.im * b.im)) / s;
	t.im = ((a.im * b.re) - (a.re * b.im)) / s;

	return t;
}

template <typename T>
__host__ __device__  struct simpleComplex<T> operator/(const struct simpleComplex<T> &a, const T &_b)
{
	struct simpleComplex<T> t, b;

    b.re = _b;
    b.im = 0.0;

	T s =  (b.re * b.re) + (b.im * b.im);

	t.re = ((a.re * b.re) + (a.im * b.im)) / s;
	t.im = ((a.im * b.re) - (a.re * b.im)) / s;

	return t;
}

__host__ __device__  float simpleComplexMod (const struct simpleComplex<float> &a)
{
	float f;

	f = sqrtf( (a.re * a.re) + (a.im * a.im) );

	return f;
}

__host__ __device__  double  simpleComplexMod (const struct simpleComplex<double> &a)
{
	double f;

	f = sqrt( (a.re * a.re) + (a.im * a.im) );

	return f;
}

template <typename T>
__host__ __device__  simpleComplex<T> sqrt(const simpleComplex<T> &a)
{
	T modval;
	simpleComplex<T> tmp;

	modval = (T)simpleComplexMod(a);

	tmp.re = sqrt( (modval + a.re) * (T)0.5 );
	tmp.im = sqrt( (modval - a.re) * (T)0.5 );

	return tmp;
}

__host__ __device__  struct simpleComplex<float> simpleComplexAdj(const struct simpleComplex<float> &a)
{
	struct simpleComplex<float> t;

	t.re = a.re;
	t.im = -a.im;

	return t;
}

__host__ __device__  struct simpleComplex<double> simpleComplexAdj(const struct simpleComplex<double> &a)
{
	struct simpleComplex<double> t;

	t.re = a.re;
	t.im = -a.im;

	return t;
}

__host__ __device__  float reciprocal(const float &a)
{
	return 1.0f/ a;
}

__host__ __device__  double reciprocal(const double &a)
{
	return 1.0 / a;
}


template <typename T>
__host__ __device__  simpleComplex<float> reciprocal(const simpleComplex<float> &a)
{
	simpleComplex<float> t;

	t.re = a.re / (a.re * a.re + a.im*a.im);
	t.im = -(a.im / (a.re * a.re + a.im*a.im));

	return t;
}

template <typename T>
__host__ __device__  simpleComplex<double> reciprocal(const simpleComplex<double> &a)
{
	simpleComplex<double> t;

	t.re = a.re / (a.re * a.re + a.im*a.im);
	t.im = -(a.im / (a.re * a.re + a.im*a.im));

	return t;
}


template <typename T>
__host__ __device__ struct simpleComplex<T> make_simpleComplex (T r, T i)
{
	simpleComplex<T> t;

	t.re = r ;
	t.im = i ;

	return t;
}
#endif


/* specialised version of float type*/

__host__ __device__ simpleComplexFloat operator+(const simpleComplexFloat &a, const simpleComplexFloat &b)
{
	simpleComplexFloat t;

	t.x = a.x + b.x;
	t.y = a.y + b.y;

	return t;
}

__host__ __device__ simpleComplexFloat operator-(const simpleComplexFloat &a, const simpleComplexFloat &b)
{
	simpleComplexFloat t;

	t.x = a.x - b.x;
	t.y = a.y - b.y;

	return t;
}

__host__ __device__ simpleComplexFloat operator+(const float &a, const simpleComplexFloat &b)
{
	simpleComplexFloat t;

	t.x = a + b.x;
	t.y = b.y;

	return t;
}

__host__ __device__  simpleComplexFloat operator-(const float &a, const simpleComplexFloat &b)
{
	simpleComplexFloat t;

	t.x = a - b.x;
	t.y = b.y;

	return t;
}

__host__ __device__  simpleComplexFloat operator+(const simpleComplexFloat &a, const float &b)
{
	simpleComplexFloat t;

	t.x = a.x + b;
	t.y = a.y ;

	return t;
}

__host__ __device__  simpleComplexFloat operator-(const simpleComplexFloat &a, const float &b)
{
	simpleComplexFloat t;

	t.x = a.x - b;
	t.y = a.y ;

	return t;
}

__host__ __device__  simpleComplexFloat operator*(const simpleComplexFloat &a, const simpleComplexFloat &b)
{
	simpleComplexFloat t;

	t.x = (a.x * b.x) - (a.y * b.y);
	t.y = (a.x * b.y) + (a.y * b.x);

	return t;
}

__host__ __device__  simpleComplexFloat operator*(const float &a, const simpleComplexFloat &b)
{
	simpleComplexFloat t;

	t.x = (a * b.x);
	t.y = (a * b.y);

	return t;
}

__host__ __device__  simpleComplexFloat operator*(const simpleComplexFloat &a, const float &b)
{
	simpleComplexFloat t;

	t.x = (a.x * b);
	t.y = (a.y * b);

	return t;
}

__host__ __device__  simpleComplexFloat operator/(const simpleComplexFloat &a, const simpleComplexFloat &b)
{
	simpleComplexFloat t;

	float s =  (b.x * b.x) + (b.y * b.y);

	t.x = ((a.x * b.x) + (a.y * b.y)) / s;
	t.y = ((a.y * b.x) - (a.x * b.y)) / s;

	return t;
}

__host__ __device__  simpleComplexFloat operator/(const float &_a, const simpleComplexFloat &b)
{
	simpleComplexFloat t, a;

	a.x = _a;
	a.y = 0.0f;

	float s =  (b.x * b.x) + (b.y * b.y);

	t.x = ((a.x * b.x) + (a.y * b.y)) / s;
	t.y = ((a.y * b.x) - (a.x * b.y)) / s;

	return t;
}

__host__ __device__  simpleComplexFloat operator/(const simpleComplexFloat &a, const float &_b)
{
	simpleComplexFloat t, b;

    b.x = _b;
    b.y = 0.0;

	float s =  (b.x * b.x) + (b.y * b.y);

	t.x = ((a.x * b.x) + (a.y * b.y)) / s;
	t.y = ((a.y * b.x) - (a.x * b.y)) / s;

	return t;
}

__host__ __device__ float simpleComplexMod (const simpleComplexFloat &a)
{
	float f;

	f = sqrtf( (a.x * a.x) + (a.y * a.y) );

	return f;
}

/*
__host__ __device__  double  simpleComplexMod (const simpleComplexDouble &a)
{
	double f;

	f = sqrt( (a.x * a.x) + (a.y * a.y) );

	return f;
}
*/

__host__ __device__  simpleComplexFloat sqrt(const simpleComplexFloat &a)
{
	float modval;
	 simpleComplexFloat tmp;

	modval = simpleComplexMod(a);

	tmp.x = sqrt( (modval + a.x) * 0.5f );
	tmp.y = sqrt( (modval - a.x) * 0.5f );

	return tmp;
}

__host__ __device__   simpleComplexFloat simpleComplexAdj(const  simpleComplexFloat &a)
{
	 simpleComplexFloat t;

	t.x = a.x;
	t.y = -a.y;

	return t;
}

/*
__host__ __device__   simpleComplexDouble simpleComplexAdj(const  simpleComplexDouble &a)
{
	 simpleComplexDouble t;

	t.x = a.x;
	t.y = -a.y;

	return t;
}
*/

__host__ __device__   simpleComplexFloat reciprocal(const  simpleComplexFloat &a)
{
	 simpleComplexFloat t;

	t.x = a.x / (a.x * a.x + a.y*a.y);
	t.y = -(a.y / (a.x * a.x + a.y*a.y));

	return t;
}

/*
__host__ __device__   simpleComplexDouble reciprocal(const  simpleComplexDouble &a)
{
	 simpleComplexDouble t;

	t.x = a.x / (a.x * a.x + a.y*a.y);
	t.y = -(a.y / (a.x * a.x + a.y*a.y));

	return t;
}
*/


__host__ __device__  simpleComplexFloat make_simpleComplexFloat (float r, float i)
{
	 simpleComplexFloat t;

	t.x = r ;
	t.y = i ;

	return t;
}


/* specialised version of double type*/

__host__ __device__ simpleComplexDouble operator+(const simpleComplexDouble &a, const simpleComplexDouble &b)
{
	simpleComplexDouble t;

	t.x = a.x + b.x;
	t.y = a.y + b.y;

	return t;
}

__host__ __device__ simpleComplexDouble operator-(const simpleComplexDouble &a, const simpleComplexDouble &b)
{
	simpleComplexDouble t;

	t.x = a.x - b.x;
	t.y = a.y - b.y;

	return t;
}

__host__ __device__ simpleComplexDouble operator+(const double &a, const simpleComplexDouble &b)
{
	simpleComplexDouble t;

	t.x = a + b.x;
	t.y = b.y;

	return t;
}

__host__ __device__  simpleComplexDouble operator-(const double &a, const simpleComplexDouble &b)
{
	simpleComplexDouble t;

	t.x = a - b.x;
	t.y = b.y;

	return t;
}

__host__ __device__  simpleComplexDouble operator+(const simpleComplexDouble &a, const double &b)
{
	simpleComplexDouble t;

	t.x = a.x + b;
	t.y = a.y ;

	return t;
}

__host__ __device__  simpleComplexDouble operator-(const simpleComplexDouble &a, const double &b)
{
	simpleComplexDouble t;

	t.x = a.x - b;
	t.y = a.y ;

	return t;
}

__host__ __device__  simpleComplexDouble operator*(const simpleComplexDouble &a, const simpleComplexDouble &b)
{
	simpleComplexDouble t;

	t.x = (a.x * b.x) - (a.y * b.y);
	t.y = (a.x * b.y) + (a.y * b.x);

	return t;
}

__host__ __device__  simpleComplexDouble operator*(const double &a, const simpleComplexDouble &b)
{
	simpleComplexDouble t;

	t.x = (a * b.x);
	t.y = (a * b.y);

	return t;
}

__host__ __device__  simpleComplexDouble operator*(const simpleComplexDouble &a, const double &b)
{
	simpleComplexDouble t;

	t.x = (a.x * b);
	t.y = (a.y * b);

	return t;
}

__host__ __device__  simpleComplexDouble operator/(const simpleComplexDouble &a, const simpleComplexDouble &b)
{
	simpleComplexDouble t;

	double s =  (b.x * b.x) + (b.y * b.y);

	t.x = ((a.x * b.x) + (a.y * b.y)) / s;
	t.y = ((a.y * b.x) - (a.x * b.y)) / s;

	return t;
}

__host__ __device__  simpleComplexDouble operator/(const double &_a, const simpleComplexDouble &b)
{
	simpleComplexDouble t, a;

	a.x = _a;
	a.y = 0.0f;

	double s =  (b.x * b.x) + (b.y * b.y);

	t.x = ((a.x * b.x) + (a.y * b.y)) / s;
	t.y = ((a.y * b.x) - (a.x * b.y)) / s;

	return t;
}

__host__ __device__  simpleComplexDouble operator/(const simpleComplexDouble &a, const double &_b)
{
	simpleComplexDouble t, b;

    b.x = _b;
    b.y = 0.0;

	double s =  (b.x * b.x) + (b.y * b.y);

	t.x = ((a.x * b.x) + (a.y * b.y)) / s;
	t.y = ((a.y * b.x) - (a.x * b.y)) / s;

	return t;
}

__host__ __device__ double simpleComplexMod (const simpleComplexDouble &a)
{
	double f;

	f = sqrt( (a.x * a.x) + (a.y * a.y) );

	return f;
}

__host__ __device__  simpleComplexDouble sqrt(const simpleComplexDouble &a)
{
	double modval;
	simpleComplexDouble tmp;

	modval = simpleComplexMod(a);

	tmp.x = sqrt( (modval + a.x) * 0.5f );
	tmp.y = sqrt( (modval - a.x) * 0.5f );

	return tmp;
}

__host__ __device__   simpleComplexDouble simpleComplexAdj(const  simpleComplexDouble &a)
{
	 simpleComplexDouble t;

	t.x = a.x;
	t.y = -a.y;

	return t;
}

/*
__host__ __device__   simpleComplexDouble simpleComplexAdj(const  simpleComplexDouble &a)
{
	 simpleComplexDouble t;

	t.x = a.x;
	t.y = -a.y;

	return t;
}
*/

__host__ __device__   simpleComplexDouble reciprocal(const  simpleComplexDouble &a)
{
	 simpleComplexDouble t;

	t.x = a.x / (a.x * a.x + a.y*a.y);
	t.y = -(a.y / (a.x * a.x + a.y*a.y));

	return t;
}

/*
__host__ __device__   simpleComplexDouble reciprocal(const  simpleComplexDouble &a)
{
	 simpleComplexDouble t;

	t.x = a.x / (a.x * a.x + a.y*a.y);
	t.y = -(a.y / (a.x * a.x + a.y*a.y));

	return t;
}
*/


__host__ __device__  simpleComplexDouble make_simpleComplexFloat (double r, double i)
{
	 simpleComplexDouble t;

	t.x = r ;
	t.y = i ;

	return t;
}

