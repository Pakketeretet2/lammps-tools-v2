#ifndef DOMAIN_HPP
#define DOMAIN_HPP


/**
   \file domain.hpp

   Declaration of domain struct and associated function.
*/


namespace lammps_tools {

/**
   \brief This struct contains information about the
          simulation box and its boundaries.
*/
struct domain {
	/// Use these bits to toggle periodicity in X, Y and Z direction.
	enum periodic_bits { BIT_X = 1,
	                     BIT_Y = 2,
	                     BIT_Z = 4 };

	double xlo[3]; ///< Lower bounds of box.
	double xhi[3]; ///< Upper bounds of box.
	int periodic;  ///< This int set periodicity, see periodic_bits

	/// Empty constructor
	domain() : xlo{0,0,0}, xhi{0,0,0}, periodic(0) {}
	/// Empty destructor
	~domain(){}

	/// Copy constructor
	domain( const domain &o );

	/// Swap, doesn't need to be friend.
	void swap( domain &f, domain &s );

	/**
	   \brief Calculates distance vector and distance^2 between two points.

	   It takes into account the periodic boundaries.

	   \param[in]   xi Point 1
	   \param[in]   xj Point 2

	   \param[out]  r  Distance vector

	   \returns     The Euclidian distance between xi and xj squared
	*/
	double dist_2( const double *xi, const double *xj, double *r ) const;

	/**
	   \brief Rewraps coordinate inside periodic boundaries:

	   \param x  Position vector to wrap
	   \param ix Image flag vector to wrap
	*/
	void rewrap( double x[3], int ix[3] ) const;

	/// Rewraps only position. \overloads rewrap.
	void rewrap( double x[3] ) const;

	/**
	   \brief Rewraps a single coordinate of a vector

	   \param x     The vector to rewrap
	   \param coord The position index to rewrap (0 for x, 1 for y, 2 for z).

	   \returns +1 if the particle moved one box image up,
	   -1 if it moved one down, 0 if it is in same box.
	*/
	template <int coord>
	int rewrap_coord( double x[3] ) const;
};

inline
double domain::dist_2( const double *xi, const double *xj, double *r ) const
{
	r[0] = xi[0] - xj[0];
	r[1] = xi[1] - xj[1];
	r[2] = xi[2] - xj[2];

	rewrap( r );

	return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}

inline
void domain::rewrap( double x[3] ) const
{
	int ix[3] = {0,0,0};
	rewrap( x, ix );
}

inline
void domain::rewrap( double x[3], int ix[3] ) const
{
	// Periodicity:
	if( periodic & BIT_X ){
		int d_ix = rewrap_coord<0>(x);
		ix[0] += d_ix;
	}
	if( periodic & BIT_Y ){
		int d_iy = rewrap_coord<1>(x);
		ix[1] += d_iy;
	}
	if( periodic & BIT_Z ){
		int d_iz = rewrap_coord<2>(x);
		ix[2] += d_iz;

	}
}

template <int coord> inline
int domain::rewrap_coord( double x[3] ) const
{
	double  L = xhi[coord] - xlo[coord];
	double  Lh = 0.5*L;
	double &rr = x[coord];

	if( rr > Lh ){
		rr -= L;
		return 1;
	}else if( rr < -Lh ){
		rr += L;
		return -1;
	}else{
		return 0;
	}
}



} // namespace lammps_tools

#endif // DOMAIN_HPP
