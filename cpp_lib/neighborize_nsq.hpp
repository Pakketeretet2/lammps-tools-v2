#ifndef NEIGHBORIZE_NSQ_HPP
#define NEIGHBORIZE_NSQ_HPP


#include <algorithm>
#include <list>
#include <vector>

#include "block_data.hpp"
#include "domain.hpp"
#include "my_assert.hpp"
#include "neighborize.hpp"

namespace lammps_tools {

namespace neighborize {

class nsq_neighborizer
{
public:
	nsq_neighborizer( const block_data &b,
	                  const std::vector<std::string> &fields,
	                  int dims, int itype, int jtype )
		: dims(dims), itype(itype), jtype(jtype), b(b),  quiet(true)

	{
		if( !grab_common_fields( b, fields, id, type, x, y, z ) ){
			my_logic_error( __FILE__, __LINE__,
			                "Failed to grab necessary fields!" );
		}
	}
	
	int dims;
	int itype, jtype;
	
	const block_data &b;
	
	std::vector<int> id;
	std::vector<int> type;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	bool quiet;

	double build( neigh_list &neighs,
	              const are_neighbours &criterion,
	              particle_filter filt = pass_all,
	              int neigh_count_estimate = 100 );
};

} // namespace neighborize

} // namespace lammps_tools


#endif // NEIGHBORIZE_NSQ_HPP
