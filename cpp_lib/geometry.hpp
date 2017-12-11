#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <algorithm>
#include <cmath>

#include "my_assert.hpp"

#include <vector>

namespace lammps_tools {

/**
   \brief Represents a point in 3D.
*/
struct point
{
	point() : x(0), y(0), z(0) {}
	point(const double *v) : x(v[0]), y(v[1]), z(v[2]) {}
	point(double x, double y, double z) : x(x), y(y), z(z) {}
	point(const std::vector<double> &v) : x(v[0]), y(v[1]), z(v[2]) {}


	point &operator=( double *v )
	{
		using std::swap;
		point n( v );
		swap( *this, n );
		return *this;
	}

	double &operator[]( int i )
	{
		my_assert( __FILE__, __LINE__, i >= 0, "Negative index not allowed!" );
		my_assert( __FILE__, __LINE__, i  < 3, "Index out of bounds!" );

		switch(i){
			case 0: return x;
			case 1: return y;
			case 2: return z;
			default: return x; // Never reached.
		}
		// Never reached:

	}
	double operator[]( int i ) const
	{
		my_assert( __FILE__, __LINE__, i >= 0, "Negative index not allowed!" );
		my_assert( __FILE__, __LINE__, i  < 3, "Index out of bounds!" );

		switch(i){
			case 0: return x;
			case 1: return y;
			case 2: return z;
		}
		// Never reached:
		return -1.337;
	}


	template <int i>
	double get() const
	{
		static_assert( i >= 0, "Negative index not allowed!" );
		static_assert( i  < 3, "Index out of bounds!" );

		switch(i){
			case 0: return x;
			case 1: return y;
			case 2: return z;
		}
		// Never reached:
		return -1.337;
	}

	double norm2() const
	{ return x*x + y*y + z*z; }

	double norm() const
	{ return sqrt(norm2()); }


	point &operator+=( const point &o )
	{
		x += o.x;
		y += o.y;
		z += o.z;
		return *this;
	}

	point &operator-=( const point &o )
	{
		x -= o.x;
		y -= o.y;
		z -= o.z;
		return *this;
	}

	point &operator*=( double f )
	{
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	point &operator/=( double f )
	{
		double invf = 1.0 / f;
		(*this) *= invf;
		return *this;
	}

	operator const double* () const
	{ return &x; }

	double x, y, z;
};

inline
point operator+( point p, point o )
{
	p += o;
	return p;
}

inline
point operator-( point p, point o )
{
	p -= o;
	return p;
}

inline
point operator*( point p, double f )
{
	p *= f;
	return p;
}

inline
point operator*( double f, point p )
{
	p *= f;
	return p;
}


inline
point operator/( point p, double f )
{
	p /= f;
	return p;
}




/**
   \brief A quaternion.
*/
struct quat {
	double a, b, c, d;

	quat() : a(0), b(0), c(0), d(0){}
	quat( double a, double b, double c, double d )
		: a(a), b(b), c(c), d(d) {}
	quat( const quat &o ) : a(o.a), b(o.b), c(o.c), d(o.d) {}
	quat( const point &p ) : a(0.0), b(p[0]), c(p[1]), d(p[2]) {}

	quat &operator=( quat o )
	{
		quat n(o);
		using std::swap;
		swap( *this, n );
		return *this;
	}

	double operator[]( int i ) const;

	template <int i>
	double get() const;

	quat &operator+=( const quat &o );
	quat &operator-=( const quat &o );
	quat &operator*=( const quat &o );

	quat &operator*=( double f );
	quat &operator/=( double f );

	quat  operator+( const quat &o ) const;
	quat  operator-( const quat &o ) const;
	quat  operator*( const quat &o ) const;
	quat  operator*( double f ) const;
	quat  operator/( double f ) const;

	quat conj() const;
	quat inv() const;


	double norm2() const;
	double norm() const;

};

inline
quat &quat::operator+=( const quat &o )
{
	a += o.a;
	b += o.b;
	c += o.c;
	d += o.d;
	return *this;
}

inline
quat &quat::operator-=( const quat &o )
{
	a -= o.a;
	b -= o.b;
	c -= o.c;
	d -= o.d;
	return *this;
}

inline
quat &quat::operator*=( const quat &o )
{
	double an = a*o.a - b*o.b - c*o.c - d*o.d;
	double bn = a*o.b + b*o.a + c*o.d - d*o.c;
	double cn = a*o.c - b*o.d + c*o.a + d*o.b;
	double dn = a*o.d + b*o.c - c*o.b + d*o.a;

	a = an;
	b = bn;
	c = cn;
	d = dn;

	return *this;
}

inline
quat quat::operator+( const quat &o ) const
{
	quat cp( a, b, c, d );
	cp += o;
	return cp;
}

inline
quat quat::operator-( const quat &o ) const
{
	quat cp( a, b, c, d );
	cp -= o;
	return cp;
}

inline
quat quat::operator*( const quat &o ) const
{
	quat cp( a, b, c, d );
	cp *= o;
	return cp;
}

inline
quat quat::conj() const
{
	return quat( a, -b, -c, -d );
}

inline
quat quat::inv() const
{
	double nn2 = norm2();
	return conj() / nn2;
}

inline
double quat::norm2() const
{
	return a*a + b*b + c*c + d*d;
}

inline
double quat::norm() const
{
	return std::sqrt( norm2() );
}

inline
double quat::operator[]( int i ) const
{
	my_assert( __FILE__, __LINE__, i >= 0, "Negative index not allowed!" );
	my_assert( __FILE__, __LINE__, i  < 4, "Index out of bounds!" );

	switch(i){
		case 0: return a;
		case 1: return b;
		case 2: return c;
		case 3: return d;
	}

	// Never reached:
	return -1.337;
}

template <int i> inline
double quat::get() const
{
	static_assert( i >= 0, "Negative index not allowed!" );
	static_assert( i  < 4, "Index out of bounds!" );

	switch(i){
		case 0: return a;
		case 1: return b;
		case 2: return c;
		case 3: return d;
	}

	// Never reached:
	return -1.337;
}

inline
quat &quat::operator*=( double f )
{
	a *= f;
	b *= f;
	c *= f;
	d *= f;
	return *this;
}

inline
quat &quat::operator/=( double f )
{
	double ff = 1.0 / f;
	a *= ff;
	b *= ff;
	c *= ff;
	d *= ff;
	return *this;
}

inline
quat  quat::operator*( double f ) const
{
	quat n(a,b,c,d);
	n *= f;
	return n;
}

inline
quat  quat::operator/( double f ) const
{
	quat n(a,b,c,d);
	n /= f;
	return n;
}

// Expand the namespace util with these dot products:
namespace util {

inline double dot( const point &p1, const point &p2 )
{
	return p1.x * p2.x + p1.y * p2.y + p1.z*p2.z;
}

inline double dot( const quat &q1, const quat &q2 )
{
	return q1.a*q2.a + q1.b*q2.b + q1.c*q2.c + q1.d*q2.d;
}

inline point cross( const point &p1, const point &p2 )
{
	return point( p1.y*p2.z - p1.z*p2.y,
	              p1.z*p2.x - p1.z*p2.z,
	              p1.x*p2.y - p1.y*p2.x );
}

}





} // lammps_tools

#endif // GEOMETRY_HPP
