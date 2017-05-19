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
   \brief Checks if special_field is set.

   \param special_field   The special field to check for.

   \returns true if the field is contained in block_data, false otherwise.
*/
bool lt_has_special_field( lt_block_data_handle *bdh, int special_field );

/**
   \brief Returns special field as vector<double>

   \param bdh            Pointer to lt_block_data_handle
   \param special_field  The special_field to get.
*/
const std::vector<double> &lt_special_field_double( lt_block_data_handle *bdh,
                                                    int special_field );


/**
   \brief Returns special field as vector<double>

   \param bdh            Pointer to lt_block_data_handle
   \param special_field  The special_field to get.
*/
const std::vector<int> &lt_special_field_int( lt_block_data_handle *bdh,
                                              int special_field );

} // extern "C"




#endif // LT_BLOCK_DATA_H
