#ifndef __cplxNum_h__
#define __cplxNum_h__

template <typename T>
struct simpleComplex {
	T re;
	T im;

	/*
	__host__ __device__  simpleComplex()
	{
		re = (T)0.0;
		im = (T)0.0;
	}

	__host__ __device__  simpleComplex(T v)
	{
		re = v;
		im = (T)0.0;
	}

	__host__ __device__  simpleComplex(T v1, T v2)
	{
		re = v1;
		im = v2;
	}
	*/
	/*
	__host__ __device__ struct simpleComplex<T>& operator += (const struct simpleComplex<T> &a)
	{
		re = re + a.re;
		im = im + a.im;

		return *this;
	}

	__host__ __device__ struct simpleComplex<T>& operator += (const T &a)
	{
		struct simpleComplex<T> t;

		re = re + a;
		//im = im;

		return *this;
	}
	*/
};


template <typename T>
__host__ __device__  struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b);

template <typename T>
__host__ __device__  struct simpleComplex<T> operator+(const T &a, const struct simpleComplex<T> &b);

template <typename T>
__host__ __device__  struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const T &b);

template <typename T>
__host__ __device__  struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b);

template <typename T>
__host__ __device__  struct simpleComplex<T> operator*(const T &a, const struct simpleComplex<T> &b);

template <typename T>
__host__ __device__  struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const T &b);

template <typename T>
__host__ __device__  struct simpleComplex<T> operator/(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b);

template <typename T>
__host__ __device__  struct simpleComplex<T> operator/(const T &_a, const struct simpleComplex<T> &b);

template <typename T>
__host__ __device__  struct simpleComplex<T> operator/(const struct simpleComplex<T> &a, const T &_b);


__host__ __device__  double simpleComplexMod (const struct simpleComplex<double> &a);
__host__ __device__  float simpleComplexMod (const struct simpleComplex<float> &a);

__host__ __device__  struct simpleComplex<float> simpleComplexAdj (const struct simpleComplex<float> &a);
__host__ __device__  struct simpleComplex<double> simpleComplexAdj (const struct simpleComplex<double> &a);


template <typename T>
__host__ __device__  T sqrt(const T &a);

template <typename T>
__host__ __device__  T reciprocal(const T &a);


template <typename T> 
__host__ __device__  struct simpleComplex<T> make_simpleComplex (T r, T i);


#endif // __cplxNum_h__
