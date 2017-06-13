#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <algorithm>
#include <cmath>

#include "my_assert.hpp"

namespace lammps_tools {

/**
   \brief Represents a point in 3D.
*/
struct point
{
	point() : x(0), y(0), z(0) {}
	explicit point(double *v) : x(v[0]), y(v[1]), z(v[2]) {}
	point(double x, double y, double z) : x(x), y(y), z(z) {}


	point &operator=( double *v )
	{
		using std::swap;
		point n( v );
		swap( *this, n );
		return *this;
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


	double x, y, z;
};

/**
   \brief A quaternion.
*/
struct quat {
	double a, b, c, d;

	quat() : a(0), b(0), c(0), d(0){}
	quat( double a, double b, double c, double d )
		: a(a), b(b), c(c), d(d) {}
	quat( const quat &o ) : a(o.a), b(o.b), c(o.c), d(o.d) {}
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







} // lammps_tools

#endif // GEOMETRY_HPP
