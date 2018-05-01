#ifndef WRITERS_HOOMD_HPP
#define WRITERS_HOOMD_HPP

/**
   \file writers_lammps.hpp

   Some writers for the HOOMD gsd format.
*/

#include <string>
#include <iosfwd>

#include <vector>

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



/**
   \brief reconstructs GSD-style buffer from vectors of data_fields.

   \param dest         Array to store the data into
   \param n_arr        The number of arrays to read in.
                       Should match field_names.size()
   \param b            The block data
   \param field_names  Names of the data fields to reconstruct


   Assumes dest is allocated and the right size.
*/
template <int data_field_type, typename T_to>
int reconstruct_fields_as_gsd( T_to *dest, int n_arr, const block_data &b,
                               const std::vector<std::string> &field_names );



} // namespace writers

} // namespace lammps_tools

#endif //WRITERS_HOOMD_HPP
