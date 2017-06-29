#include "msd.hpp"

namespace lammps_tools {

namespace fluctuations {


double compute_msd( readers::dump_reader *reader, std::vector<double> &msd,
                    std::list<bigint> &time, const std::list<int> &ids,
                    bool got_ids )
{
	return 0.0;
}

} // namespace fluctuations

} // namespace lammps_tools
