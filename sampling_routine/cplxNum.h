#ifndef __cplxNum_h__
#define __cplxNum_h__

#include <cuComplex.h>

#if 0
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
#endif

typedef cuFloatComplex simpleComplexFloat;

/*
struct  __align__(8) simpleComplexFloat {
	float re;
	float im;
};

struct __align__(8) simpleComplexDouble {
	double re;
	double im;
};
*/

struct __align__(16) simpleComplexDouble {
	double x;
	double y;
};


#if 0

//typedef struct __simpleComplexFloat simpleComplexFloat;
//typedef struct __simpleComplexDouble simpleComplexDouble;

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

#endif

/* specialised version of float type*/

__host__ __device__  simpleComplexFloat operator+(const simpleComplexFloat &a, const simpleComplexFloat &b);

__host__ __device__  simpleComplexFloat operator+(const float &a, const  simpleComplexFloat &b);

__host__ __device__  simpleComplexFloat operator+(const  simpleComplexFloat &a, const float &b);

__host__ __device__  simpleComplexFloat operator*(const  simpleComplexFloat &a, const  simpleComplexFloat &b);

__host__ __device__  simpleComplexFloat operator*(const float &a, const  simpleComplexFloat &b);

__host__ __device__  simpleComplexFloat operator*(const  simpleComplexFloat &a, const float &b);

__host__ __device__  simpleComplexFloat operator/(const  simpleComplexFloat &a, const  simpleComplexFloat &b);

__host__ __device__  simpleComplexFloat operator/(const float &_a, const  simpleComplexFloat &b);

__host__ __device__  simpleComplexFloat operator/(const  simpleComplexFloat &a, const float &_b);


__host__ __device__  float simpleComplexMod (const  simpleComplexFloat &a);
//__host__ __device__  double simpleComplexMod (const  simpleComplexDouble &a);

__host__ __device__   simpleComplexFloat simpleComplexAdj (const  simpleComplexFloat &a);
//__host__ __device__   simpleComplexDouble simpleComplexAdj (const  simpleComplexDouble &a);


__host__ __device__   simpleComplexFloat sqrt(const  simpleComplexFloat &a);
//__host__ __device__   simpleComplexDouble reciprocal(const  simpleComplexDouble &a);

//__host__ __device__  float sqrt(const float &a);
//__host__ __device__  float reciprocal(const float &a);

//__host__ __device__  double sqrt(const double &a);
//__host__ __device__  double reciprocal(const double &a);

__host__ __device__   simpleComplexFloat make_simpleComplexFloat (float r, float i);
//__host__ __device__   simpleComplexDouble make_simpleComplexDouble (double r, double i);


/* specialised version of double type*/

__host__ __device__  simpleComplexDouble operator+(const simpleComplexDouble &a, const simpleComplexDouble &b);

__host__ __device__  simpleComplexDouble operator+(const double &a, const  simpleComplexDouble &b);

__host__ __device__  simpleComplexDouble operator+(const  simpleComplexDouble &a, const double &b);

__host__ __device__  simpleComplexDouble operator*(const  simpleComplexDouble &a, const  simpleComplexDouble &b);

__host__ __device__  simpleComplexDouble operator*(const double &a, const  simpleComplexDouble &b);

__host__ __device__  simpleComplexDouble operator*(const  simpleComplexDouble &a, const double &b);

__host__ __device__  simpleComplexDouble operator/(const  simpleComplexDouble &a, const  simpleComplexDouble &b);

__host__ __device__  simpleComplexDouble operator/(const double &_a, const  simpleComplexDouble &b);

__host__ __device__  simpleComplexDouble operator/(const  simpleComplexDouble &a, const double &_b);


__host__ __device__  double simpleComplexMod (const  simpleComplexDouble &a);
//__host__ __device__  double simpleComplexMod (const  simpleComplexDouble &a);

__host__ __device__   simpleComplexDouble simpleComplexAdj (const  simpleComplexDouble &a);
//__host__ __device__   simpleComplexDouble simpleComplexAdj (const  simpleComplexDouble &a);


__host__ __device__   simpleComplexDouble sqrt(const  simpleComplexDouble &a);
//__host__ __device__   simpleComplexDouble reciprocal(const  simpleComplexDouble &a);

//__host__ __device__  float sqrt(const float &a);
//__host__ __device__  float reciprocal(const float &a);

//__host__ __device__  double sqrt(const double &a);
//__host__ __device__  double reciprocal(const double &a);

__host__ __device__   simpleComplexDouble make_simpleComplexDouble (double r, double i);


#endif // __cplxNum_h__
