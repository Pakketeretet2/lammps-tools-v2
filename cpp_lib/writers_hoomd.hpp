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

enum gsd_write_props {
	TIME_STEP       = 1 << 0,
	DIMENSIONS      = 1 << 1,
	BOX             = 1 << 2,
	PARTICLE_NUMBER = 1 << 3,
	POSITIONS       = 1 << 4,
	TYPEID          = 1 << 5,
	TYPES           = 1 << 6,
	ORIENTATION     = 1 << 7,
	BODY            = 1 << 8,
	MOM_INERTIA     = 1 << 9,

	ALL_PROPS       = (1 << 10) - 1 // Dummy to encode for all.
};

int block_to_hoomd_gsd( const std::string &fname, const block_data &b,
                        const std::string &write_mode,
                        uint props = ALL_PROPS );

/**
   \brief Write block data to a gsd file through a gsd_handle.

   \param gh     Ptr to the handle to the gsd file
   \param b      Block data to write
   \param props  A bitfield that specifies which things to write to file

   \returns a status (0 on success, non-zero otherwise).
*/
int block_to_hoomd_gsd( gsd_handle *gh, const block_data &b,
                        uint props = ALL_PROPS );

} // namespace writers

} // namespace lammps_tools

#endif //WRITERS_HOOMD_HPP
