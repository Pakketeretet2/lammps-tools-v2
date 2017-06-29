#ifndef NEIGHBORIZE_NSQ_HPP
#define NEIGHBORIZE_NSQ_HPP


#include <algorithm>
#include <list>
#include <vector>

#include "block_data.hpp"
#include "domain.hpp"
#include "my_assert.hpp"
#include "neighborize.hpp"

/**
   \file neighborize_nsq.hpp
*/


namespace lammps_tools {

namespace neighborize {

class neighborizer_nsq : public neighborizer
{
public:
	neighborizer_nsq( const block_data &b, const std::vector<int> &s1,
	                  const std::vector<int> &s2, int dims )
		: neighborizer( b, s1, s2, dims )
	{}

	~neighborizer_nsq(){}

private:
	virtual int build( neigh_list &neighs,
	                   const are_neighbours &criterion );
};

} // namespace neighborize

} // namespace lammps_tools


#endif // NEIGHBORIZE_NSQ_HPP
