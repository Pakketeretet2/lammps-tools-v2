#ifndef LT_BLOCK_DATA_H
#define LT_BLOCK_DATA_H

#include "../cpp_lib/my_assert.hpp"
#include "../cpp_lib/types.hpp"


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
	lt_block_data_handle( const lt_block_data_handle &o )
	{
		bd = new lammps_tools::block_data;
		*bd = *o.bd;
	}

	lt_block_data_handle &operator=( lt_block_data_handle o )
	{
		using std::swap;
		lt_block_data_handle new_block( o );
		swap( *this, new_block );
		return *this;
	}
	
	~lt_block_data_handle()
	{
		std::cerr << "Calling delete on contained block_data at "
		          << bd << ".\n";
		if( bd ) delete bd;
	}

	lammps_tools::bigint time_step() const
	{ return bd->tstep; }
	
	lammps_tools::bigint n_atoms() const
	{ return bd->N; }

	int n_atom_types() const
	{ return bd->tstep; }

	int atom_style() const
	{ return bd->atom_style; }

	lammps_tools::block_data *bd;
};


/**
   \brief returns a pointer to named data field.
   
   \param[in]  bdh  The block data to read from.
   \param[in]  name The name of the data to read.
   \param[out] size Will contain the data field size on success.
   \param[out] type Will contain the data type on success.

   On failure to find data, size and type shall be unchanged.

   \returns a pointer to named data field or nullptr if it could not be found.
*/
lammps_tools::data_field *lt_get_data( lt_block_data_handle *bdh,
                                       const char *name,
                                       size_t *size, int *type );

/**
   \brief returns a pointer to named data field, interpreted as double.
   
   \param[in]  bdh  The block data to read from.
   \param[in]  name The name of the data to read.
   \param[out] size Will contain the data field size on success.
   \param[out] type Will contain the data type on success.

   On failure to find data, size, type and data shall be unchanged.

   \returns true on success, false otherwise.
*/
bool lt_data_as_double( lt_block_data_handle *bdh, const char *name,
                        size_t *size, int *type, std::vector<double> &data );

/*
bool lt_data_as_double( lt_block_data_handle *bdh, const char *name,
                        size_t *size, int *type, double **data )
*/

/**
   \brief returns a pointer to named data field, interpreted as int.
   
   \param[in]  bdh  The block data to read from.
   \param[in]  name The name of the data to read.
   \param[out] size Will contain the data field size on success.
   \param[out] type Will contain the data type on success.

   On failure to find data, size and type shall be unchanged.

   \returns a pointer to named data field or nullptr if it could not be found.
*/
/*
lammps_tools::data_field_int *lt_data_as_int( lt_block_data_handle *bdh,
                                              const char *name,
                                              size_t *size, int *type );
*/


/**
   \brief returns an empty block_data_handle.

   \returns an empty block_data_handle.
*/
lt_block_data_handle lt_empty_block_data_handle();

} // extern "C"




#endif // LT_BLOCK_DATA_H
