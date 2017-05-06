#ifndef LT_DUMP_READER_H
#define LT_DUMP_READER_H

/**
   \file lt_dump_reader.h

   This part of the interface exposes a dump reader via a handle.
*/

#include "../cpp_lib/dump_reader.hpp"

#include "lt_block_data.h"

extern "C" {

/**
   \brief Enumerates possible dump reader statuses.
*/
enum LT_DUMP_READER_STATUS {
	IS_GOOD = 0,
	AT_EOF  = 1,
	IS_BAD  = -1,
	POINTER_NULL = -2
};

/**
   \brief contains a pointer to an instantiated dump reader.
*/
struct lt_dump_reader_handle
{
	lammps_tools::readers::dump_reader *dr;
	int fformat, dformat;
};

// Expose the interface of dump_reader through functions and handles:


/**
   \brief Creates a new dump_reader and returns the handle.

   For the recognised file and dump fomats, see FILE_FORMATS and
   DUMP_FORMATS in cpp_lib/dump_reader.hpp.

   \param fname    Dump file name to open. Use "-" for stdin.
   \param fformat  Specifies the file format.
   \param dformat  Specifies the dump format.

   \returns a handle to a new dump reader handle. The internal
            dump_reader pointer can be nullptr. Use
*/
lt_dump_reader_handle lt_new_dump_reader( const char *fname,
                                          int fformat, int dformat );

/**
   \brief Cleans up and releases given dump reader handle.

   \param drh   Dump reader handle to release.
*/
void lt_delete_dump_reader( lt_dump_reader_handle drh );

/**
   \brief Returns the status of the dump reader handle.

   See lt_dump_reader_status for more info.

   \param drh The dump reader handle to check.

   \returns an int describing the status of drh.
*/
int lt_dump_reader_status( lt_dump_reader_handle drh );


/**
   \brief Reads the next block_data object from given dump_reader.

   \param[in]   drh The lt_dump_reader_handle to check.
   \param[out]  bdh The lt_block_data_handle to write the block data to.

   \returns An int describing the status of drh, see lt_dump_reader_status.
*/
int lt_get_next_block( lt_dump_reader_handle drh, lt_block_data_handle *bdh );


/**
   \brief Counts the number of blocks left in the given dump file.

   \param[in/out] drh The lt_dump_reader_handle to check.

   \warning       The state dump_reader_handle is modified by this function!
                  If you just want to know how much blocks a given dump file
                  has, create a fresh dump_reader, call this function, then
                  create another fresh dump_reader to do analysis with.

   \returns An int describing the status of drh, see lt_dump_reader_status.
*/
int lt_number_of_blocks( lt_dump_reader_handle drh );


/**
   \brief Sets the column headers that some dump formats expect for dump reader.

   \param drh      Handle to the dump reader.
   \param n        The column index to set.
   \param headers  The header as array of '\0'-terminated char arrays.
*/
void lt_set_col_header( lt_dump_reader_handle drh, int n, const char *header );


} // extern "C"

#endif // LT_DUMP_READER_H