#ifndef RMSD_FLUCTUATIONS_HPP
#define RMSD_FLUCTUATIONS_HPP

/**
   \file rmsd_fluctuations.hpp

   Contains some routines for calculating root-mean-square fluctuations,
   normal mode analysis, etc.
*/

#include <array>
#include <list>
#include <vector>

#include "dump_reader.hpp"
#include "types.hpp"

namespace lammps_tools {

namespace fluctuations {

using lammps_tools::bigint;
typedef std::array<double,4> vec4;

/**
   \brief Calculates the rood mean square deviation of atoms in given dump file
          with respect to a running average of their positions.

   \param[in/out] reader  The dump reader whose time steps to read
   \param[out]    rmsd    Will contain the root-mean-square fluctuations
   \param[out]    time    Will contain the time steps corresponding to rmsd
   \param[in]     ids     Only calculate the RMSD for these ids.
*/
double rmsd( readers::dump_reader *reader, std::list<double> &rmsds,
             std::list<bigint> &time, const std::list<int> &ids,
             bool got_ids = true );



/**
   \brief Calculates the rood mean square deviation of atoms in given dump file
          with respect to a running average of their positions.

   \param[in/out] reader  The dump reader whose time steps to read
   \param[out]    rmsd    Will contain the root-mean-square fluctuations
   \param[out]    time    Will contain the time steps corresponding to rmsd
*/
double rmsd( readers::dump_reader *reader, std::list<double> &rmsds,
             std::list<bigint> &time );




} // namespace fluctuations

} // namespace lammps_tools




#endif // RMSD_FLUCTUATIONS_HPP
