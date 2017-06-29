#ifndef MSD_HPP
#define MSD_HPP


#include "readers.hpp"
#include "types.hpp"

#include <list>
#include <vector>

namespace lammps_tools {

namespace fluctuations {


double compute_msd( readers::dump_reader *reader, std::vector<double> &msd,
                    std::list<bigint> &time, const std::list<int> &ids,
                    bool got_ids );

} // namespace fluctuations

} // namespace lammps_tools


#endif // MSD_HPP
