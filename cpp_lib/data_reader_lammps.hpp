#ifndef DATA_READER_LAMMPS_HPP
#define DATA_READER_LAMMPS_HPP

#include <iosfwd>
#include <string>

#include "block_data.hpp"

/// Contains functions and classes that are related to reading data files.
namespace data_readers {

/**
   \brief Reads block from input stream.

   \param in      The input file stream.
   \param status  Set to 0 on success, non-negative otherwise.
   
   \returns       A new block_data object.
*/
block_data block_data_from_lammps_data( std::istream &in, int &status,
                                        bool quiet = true );


} // namespace data_readers

#endif // DATA_READER_LAMMPS_HPP
