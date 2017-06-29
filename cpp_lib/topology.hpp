#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

/**
   \file topology.hpp

   Declarations/definitions related to molecular topology.
*/

#include "types.hpp"

#include <vector>

namespace lammps_tools {

// Need to forward-declare:
class block_data;

/**
   \brief Defines a general molecular topology.
*/
struct topology
{
	/// Swaps \p f and \p s.
	friend void swap( topology &f, topology &s );
};

/**
   \brief Defines a bond between atoms.
*/
struct bond {
	bigint id;   ///< ID of bond.
	int type; ///< Type of bond.
	int particle1, particle2; ///< Particles in bond, by index.
};

/**
   \brief Defines an angle between atom pairs.
*/
struct angle {
	bigint id;   ///< ID of bond.
	int type; ///< Type of angle.
	int particle1, particle2, particle3; ///< Particles in angle, by index.
};



}

#endif // TOPOLOGY_HPP
