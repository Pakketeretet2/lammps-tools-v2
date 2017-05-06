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
	
};

inline
double domain::dist_2( const double *xi, const double *xj, double *r ) const
{
	r[0] = xi[0] - xj[0];
	r[1] = xi[1] - xj[1];
	r[2] = xi[2] - xj[2];

	double L[3];
	double Lh[3];
	L[0] = xhi[0] - xlo[0];
	L[1] = xhi[1] - xlo[1];
	L[2] = xhi[2] - xlo[2];

	Lh[0] = 0.5*L[0];
	Lh[1] = 0.5*L[1];
	Lh[2] = 0.5*L[2];

	auto rewrap = [&r, L, Lh]( int d ) -> void {
		if( r[d] > Lh[d] ){
			r[d] -= L[d];
		}else if( r[d] < -Lh[d] ){
			r[d] += L[d];
		} };
	
	// Periodicity:
	if( periodic & BIT_X ){
		rewrap(0);
	}
	if( periodic & BIT_Y ){
		rewrap(1);
	}
	if( periodic & BIT_Z ){
		rewrap(2);
	}

	return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}
	



} // namespace lammps_tools

#endif // DOMAIN_HPP
