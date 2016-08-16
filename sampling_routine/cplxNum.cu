#include "cplxNum.cuh"

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
