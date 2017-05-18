#ifndef WRITERS_LAMMPS_HPP
#define WRITERS_LAMMPS_HPP

/**
   \file writers_lammps.hpp

   Some writers for lammps format.
*/

#include <string>
#include <iosfwd>
#include "block_data.hpp"
#include "dump_reader.hpp"

namespace lammps_tools {

/// \brief Contains functions for writing data.
namespace writers {

/**
   \brief Writes block_data to LAMMPS data file.

   \param fname  output file name
   \param b      block_data to write
*/
void block_to_lammps_data( const std::string &fname, const block_data &b );

/**
   \brief Writes block_data to output stream in LAMMPS data format.

   \overloads block_to_lammps_data
*/
void block_to_lammps_data( std::ostream &out, const block_data &b );


/**
   \brief Writes block_data to LAMMPS dump file.

   \param fname    output file name
   \param b        block_data to write
   \param fformat  file format to use (see lammps_tools::readers::FILE_FORMATS
                   in dump_readers.hpp)
*/
void block_to_lammps_dump( const std::string &fname, const block_data &b,
                           int fformat );

/**
   \brief Writes block_data to output stream in LAMMPS dump format.

   \overloads block_to_lammps_dump
*/
void block_to_lammps_dump( std::ostream &out, const block_data &b, int fformat );

/// Writes block to output stream in plain text LAMMPS dump format.
void block_to_lammps_dump_text( std::ostream &out, const block_data &b );

/// Writes block to output stream in binary LAMMPS dump format.
void block_to_lammps_dump_bin( std::ostream &out, const block_data &b );


} // namespace writers

} // namespace lammps_tools

#endif // WRITERS_LAMMPS_HPP
