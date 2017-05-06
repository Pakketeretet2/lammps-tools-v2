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

namespace writers {

void block_to_lammps_data( const std::string &fname, const block_data &b );
void block_to_lammps_data( std::ostream &out, const block_data &b );

void block_to_lammps_dump( const std::string &fname, const block_data &b,
                           int fformat );

void block_to_lammps_dump( std::ostream &out, const block_data &b, int fformat );

void block_to_lammps_dump_text( std::ostream &out, const block_data &b );

void block_to_lammps_dump_bin( std::ostream &out, const block_data &b );


} // namespace writers

} // namespace lammps_tools

#endif // WRITERS_LAMMPS_HPP
