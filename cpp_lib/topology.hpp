#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

/**
   \file topology.hpp
   
   Declarations/definitions related to molecular topology.
*/

namespace lammps_tools {

struct topology
{
	/// Swaps \p f and \p s.
	friend void swap( topology &f, topology &s );
};

}

#endif // TOPOLOGY_HPP
