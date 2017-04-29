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
	domain(){}
	/// Empty destructor
	~domain(){}

	/// Copy constructor
	domain( const domain &o );

	/// Swap, doesn't need to be friend. 
	void swap( domain &f, domain &s );
	
};

} // namespace lammps_tools

#endif // DOMAIN_HPP
