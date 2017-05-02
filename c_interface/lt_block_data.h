#ifndef LT_BLOCK_DATA_H
#define LT_BLOCK_DATA_H

#include "../cpp_lib/my_assert.hpp"

/**
   \file lt_block_data.h
   
   Exposes a block_data object through the C interface.
*/

#include "../cpp_lib/block_data.hpp"

extern "C" {

/**
   \brief contains a pointer to an instantiated dump reader.
*/
struct lt_block_data_handle
{
	lt_block_data_handle() : bd(nullptr)
	{
		bd = new lammps_tools::block_data;
	}
	~lt_block_data_handle()
	{
		if( bd ) delete bd;
	}
	
	lammps_tools::block_data *bd;
};


/**
   \brief extracts the named data field from given block data object.

   \param[in]    bdf   The block data object whose data to access.
   \param[in]    name  Name of the data field.
   \param[out]   size  The size of the array pointed to.

   \returns a pointer to the data, or nullptr if the data does not exist.
*/
void *lt_get_data( lt_block_data_handle bdh, const char *name, size_t *size );


} // extern "C"




#endif // LT_BLOCK_DATA_H
