#ifndef LT_BLOCK_WRITERS_H
#define LT_BLOCK_WRITERS_H

#include "lt_block_data.h"
#include "../cpp_lib/writers.hpp"

#include <iosfwd>

extern "C" {

/**
   \brief Writes block as lammps data to given file name.

   \param[in] fname  The file name to write
   \param[in] w_mode String that contains the file write mode
   \param[in] bdh    block_data_handle to write
*/
int lt_block_writers_lammps_data( const char *fname, const char *w_mode,
                                  const lt_block_data_handle *bdh );

int lt_block_writers_lammps_dump( const char *fname, const char *w_mode,
                                  const lt_block_data_handle *bdh );

int lt_block_writers_hoomd_gsd( const char *fname, const char *w_mode,
                                const lt_block_data_handle *bdh );


} // extern "C"

int c_like_open_out_file( const char *fname, const char *w_mode,
                          std::ofstream &out, std::ios_base::openmode &wmode );

#endif // LT_BLOCK_WRITERS_H
