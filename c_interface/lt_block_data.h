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
		/*
		std::cerr << "This is the destructor of " << this
		          << ", I am going to delete data @" << bd << "\n";
		*/
		if( bd ) delete bd;
	}

	lammps_tools::bigint time_step() const
	{ return bd->tstep; }

	lammps_tools::bigint n_atoms() const
	{ return bd->N; }

	void set_n_atoms( lammps_tools::bigint N )
	{ bd->set_natoms( N ); }

	void set_atom_style( int atom_style )
	{ bd->atom_style = atom_style; }

	int n_atom_types() const
	{ return bd->tstep; }

	int atom_style() const
	{ return bd->atom_style; }

	int n_types() const
	{ return bd->N_types; }

	void set_n_types( int N )
	{ bd->set_ntypes( N ); }

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

/**
   \brief Sets the meta-data for the block_data.

*/
void lt_block_data_set_meta( lt_block_data_handle *bdh,
                             lammps_tools::bigint tstep,
                             lammps_tools::bigint natoms,
                             double xlo, double xhi,
                             double ylo, double yhi,
                             double zlo, double zhi,
                             int periodic_bits, int atom_style );

void lt_block_data_set_domain( lt_block_data_handle *bdh,
                               double xlo, double xhi,
                               double ylo, double yhi,
                               double zlo, double zhi,
                               int periodic_bits );


void lt_block_data_print_stats( lt_block_data_handle *bdh );

void lt_block_data_remove_field( lt_block_data_handle *bdh, const char *name );
void lt_block_data_swap_fields( lt_block_data_handle *bdh, const char *name,
                                const lt_data_field_handle *df );
int lt_block_data_get_domain_periodic( const lt_block_data_handle *bdh );
const double *lt_block_data_get_domain_xlo( const lt_block_data_handle *bdh );
const double *lt_block_data_get_domain_xhi( const lt_block_data_handle *bdh );

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

/**
   \brief Creates a new, empty block_data_handle.
*/
lt_block_data_handle *lt_new_block_data_handle();

/**
   \brief Deletes a block_data_handle
*/
void lt_delete_block_data_handle( lt_block_data_handle *bdh );


/**
   \brief Stores filtered block_data in given handle.

   \warning This function deletes whatever was in the old handle!
*/
void lt_block_data_filter( lt_block_data_handle *dest, int size, const void *ids,
                           const lt_block_data_handle *src );



const std::vector<double>
lt_block_data_get_domain_xlo_vec( const lt_block_data_handle *bdh );

const std::vector<double>
lt_block_data_get_domain_xhi_vec( const lt_block_data_handle *bdh );




#endif // LT_BLOCK_DATA_H
