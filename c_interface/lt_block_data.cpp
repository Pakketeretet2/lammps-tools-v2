#include "lt_block_data.h"


lt_block_data_handle lt_empty_block_data_handle()
{
	lt_block_data_handle b;
	return b;
}

bool lt_has_special_field( lt_block_data_handle *bdh, int special_field )
{
	return bdh->bd->get_special_field( special_field ) != nullptr;
}



const std::vector<double> &lt_special_field_double( lt_block_data_handle *bdh,
                                                    int special_field )
{
	const lammps_tools::data_field *df
		= bdh->bd->get_special_field( special_field );
	if( !df ){
		std::cerr << "No data field for special field type "
		          << special_field << "!\n";
		std::terminate();
	}
	return lammps_tools::data_as<double>( df );
}


const std::vector<int> &lt_special_field_int( lt_block_data_handle *bdh,
                                              int special_field )
{
	const lammps_tools::data_field *df
		= bdh->bd->get_special_field( special_field );
	if( !df ){
		std::cerr << "No data field for special field type "
		          << special_field << "!\n";
		std::terminate();
	}
	return lammps_tools::data_as<int>( df );
}


lt_data_field_handle lt_data_by_name( lt_block_data_handle *bdh, const char *name )
{
	const std::string n(name);
	int i = 0;

	lt_data_field_handle ldf;
	ldf.df = nullptr;
	const lammps_tools::block_data &b = *bdh->bd;

	for( int i = 0; i < b.n_data_fields(); ++i ){
		if( n == b[i].name ){
			const lammps_tools::data_field *df = &b[i];
			ldf.df = df;
			break;
		}
	}
	return ldf;
}


lt_data_field_handle lt_data_by_index( lt_block_data_handle *bdh, int i )
{
	lt_data_field_handle ldf;
	ldf.df = nullptr;
	const lammps_tools::block_data &b = *bdh->bd;

	if( i < 0 || i >= b.n_data_fields() ){
		return ldf;
	}
	ldf.df = &b[i];
	return ldf;
}


int lt_n_data_fields( lt_block_data_handle *bdh )
{
	return bdh->bd->n_data_fields();
}
