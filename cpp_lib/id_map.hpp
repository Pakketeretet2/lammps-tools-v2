#ifndef ID_MAP_HPP
#define ID_MAP_HPP

/**
   \file id_map.hpp
   
   Declarations for a mapping of atom ids to atom indices.
*/

#include <algorithm>
#include <vector>

#include "my_assert.hpp"

namespace lammps_tools {

/**
   Contains a mapping from atom ids to indices.
*/
class id_map {
public:
	/**
	   Constructor based on given atom ids

	   \param ids  Vector containing the atom ids.
	*/
	template <typename int_type> explicit
	id_map( const std::vector<int_type> &ids )
	{
		auto max_id_p = std::max_element( ids.begin(), ids.end() );
		my_assert( __FILE__, __LINE__, max_id_p != ids.end(),
		           "Ids was empty!" );
		int_type max_id = *max_id_p;
		m.resize( max_id + 1 );
		for( int i = 0; i < ids.size(); ++i ){
			m[ids[i]] = i;
		}

	}

	/// Returns the index in ids (and other arrays) of given id
	int operator[]( std::size_t id ) const
	{
		return m[id];
	}
	
private:
	std::vector<int> m; ///< std::map the mapping is stored in
};

} // namespace lammps_tools


#endif // ID_MAP_HPP
