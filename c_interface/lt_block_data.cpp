#include "lt_block_data.h"
// #include "../cpp_lib/block_data_access.hpp"

lt_block_data_handle lt_empty_block_data_handle()
{
	// Cannot throw:
	lt_block_data_handle b;
	return b;
}

bool lt_has_special_field( lt_block_data_handle *bdh, int special_field )
{
	// Cannot throw:
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

	lt_data_field_handle dfh;
	dfh.df = df;
	return lt_data_as_double_vec( dfh );
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
	lt_data_field_handle dfh;
	dfh.df = df;
	return lt_data_as_int_vec( dfh );
}


lt_data_field_handle lt_data_by_name( lt_block_data_handle *bdh, const char *name )
{
	const std::string n(name);
	lt_data_field_handle ldf;
	ldf.df = nullptr;
	const lammps_tools::block_data &b = *bdh->bd;

	for( std::size_t i = 0; i < b.n_data_fields(); ++i ){
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
	std::size_t ii = i;
	if( i < 0 || ii >= b.n_data_fields() ){
		return ldf;
	}else{
		std::cerr << "Index " << i << " is out of range! "
		          << "Ignoring call to lt_data_by_index!\n";
	}
	ldf.df = &b[i];
	return ldf;
}


int lt_n_data_fields( lt_block_data_handle *bdh )
{
	return bdh->bd->n_data_fields();
}
