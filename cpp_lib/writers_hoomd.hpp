#ifndef WRITERS_HOOMD_HPP
#define WRITERS_HOOMD_HPP

/**
   \file writers_lammps.hpp

   Some writers for the HOOMD gsd format.
*/

#include <string>
#include <iosfwd>

struct gsd_handle;

namespace lammps_tools {

class block_data;

/// \brief Contains functions for writing data.
namespace writers {

int block_to_hoomd_gsd( const std::string &fname, const block_data &b,
                        const std::string &write_mode = "wb" );

int block_to_hoomd_gsd( gsd_handle *gh, const block_data &b );

} // namespace writers

} // namespace lammps_tools

#endif //WRITERS_HOOMD_HPP
