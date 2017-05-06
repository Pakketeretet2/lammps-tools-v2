#ifndef DATA_READER_LAMMPS_HPP
#define DATA_READER_LAMMPS_HPP

/**
   \file data_reader_lammps.hpp

   \brief Contanins functions for reading in data files.
*/

#include <iosfwd>
#include <string>

#include "block_data.hpp"

namespace lammps_tools {

/// Contains functions and classes that are related to reading data files.
namespace readers {

/**
   \brief Reads block from input stream.

   \param in      The input file stream.
   \param status  Set to 0 on success, non-negative otherwise.

   \returns       A new block_data object.
*/
block_data block_data_from_lammps_data( std::istream &in, int &status,
                                        bool quiet = true );

} // namespace readers

} // namespace lammps_tools

#endif // DATA_READER_LAMMPS_HPP
