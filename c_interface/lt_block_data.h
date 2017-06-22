#ifndef LT_BLOCK_DATA_H
#define LT_BLOCK_DATA_H

#include "../cpp_lib/my_assert.hpp"
#include "../cpp_lib/types.hpp"

#include "lt_data_field.h"

/**
   \file lt_block_data.h

   Exposes a block_data object through the C interface.
*/

#include "../cpp_lib/block_data.hpp"
#include "../cpp_lib/util.hpp"

#include <memory>

extern "C" {

/**
   \brief contains a pointer to a block_data, either empty or initialised.
*/
struct lt_block_data_handle
{
	lt_block_data_handle() : bd(nullptr)
	{
		using lammps_tools::util::make_unique;
		bd = new lammps_tools::block_data;
	}

	lt_block_data_handle( const lt_block_data_handle &o ) : bd(nullptr)
	{
		bd = new lammps_tools::block_data;
		// No worries, this calls the copy constructor of block_data:
		*bd = *o.bd;
	}

	lt_block_data_handle &operator=( const lt_block_data_handle &o )
	{
		// I think this might break now, you need to copy or move.
		using std::swap;
		lt_block_data_handle new_block( o );
		swap( *this, new_block );
		// bd->swap( new_block.bd );
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

	int n_types() const
	{ return bd->N_types; }

	const lammps_tools::block_data &get_const_ref() const
	{ return *bd; }

	lammps_tools::block_data *get_ptr()
	{ return bd; }

	lammps_tools::block_data *bd;
};


/**
   \brief Checks if special_field is set.

   \param special_field   The special field to check for.

   \returns true if the field is contained in block_data, false otherwise.
*/
bool lt_has_special_field( lt_block_data_handle *bdh, int special_field );


/**
   \brief Returns data field specified by name, or nullptr if it doesn't exist.
*/
int lt_n_data_fields( lt_block_data_handle *bdh );


/**
   \brief Returns data field specified by name, or nullptr if it doesn't exist.
*/
lt_data_field_handle lt_data_by_name( lt_block_data_handle *bdh, const char *name );

/**
   \brief Returns data field by index, or nullptr if index out of range
*/
lt_data_field_handle lt_data_by_index( lt_block_data_handle *bdh, int i );


/**
   \brief Adds given data field to block_data.

   \param[in/out] bdh  Block_data to add data to.
   \param[in]     dfh  Data field to add. This is _copied_.
*/
void lt_block_data_add_data_field( lt_block_data_handle *bdh,
                                   const lt_data_field_handle *dfh );


/**
   \brief Adds given data field to block_data as a special field.

   \param[in/out] bdh  Block_data to add data to.
   \param[in]     dfh  Data field to add. This is _copied_.
   \parma[in]     type The type of the special field.
*/
void lt_block_data_add_special_field( lt_block_data_handle *bdh,
                                      const lt_data_field_handle *dfh,
                                      int type );

} // extern "C"

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





#endif // LT_BLOCK_DATA_H
