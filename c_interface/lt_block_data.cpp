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
