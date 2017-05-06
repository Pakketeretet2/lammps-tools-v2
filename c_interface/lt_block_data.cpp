#include "lt_block_data.h"


lammps_tools::data_field *lt_get_data( lt_block_data_handle *bdh,
                                       const char *name,
                                       size_t *size, int *type )
{
	lammps_tools::data_field *df = bdh->bd->get_data_rw( name );
	if( !df ){
		return nullptr;
	}
	*type = df->type();
	*size = df->size();
	
	return df;
}


lt_block_data_handle lt_empty_block_data_handle()
{
	lt_block_data_handle b;
	return b;
}

bool lt_data_as_double( lt_block_data_handle *bdh, const char *name,
                        size_t *size, int *type, double **data )
{
	return false;
}

bool lt_data_as_double( lt_block_data_handle *bdh, const char *name,
                        size_t *size, int *type, std::vector<double> &data )
{
	int old_size = *size;
	lammps_tools::data_field *df = lt_get_data( bdh, name, size, type );
	if( !df ) return false;

	if( *type != lammps_tools::data_field::DOUBLE ){
		*size = old_size;
		return false;
	}

	using dfd = lammps_tools::data_field_double;
	dfd *d = static_cast<dfd*>( df );
	data = d->get_data_rw();
	
	return true;
}
