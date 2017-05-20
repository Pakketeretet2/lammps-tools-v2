#include "lt_data_field.h"

#include "../cpp_lib/data_field.hpp"

int lt_data_field_size( lt_data_field_handle d )
{
	return d.df->size();
}

const double *lt_data_as_double( lt_data_field_handle d )
{
	using dfd = lammps_tools::data_field_double;
	if( d.df->type() != lammps_tools::data_field::DOUBLE ){
		return nullptr;
	}else{
		const dfd* dd = static_cast<const dfd*>( d.df );
		return dd->get_data_ptr();
	}
}

const int *lt_data_as_int( lt_data_field_handle d )
{
	using dfi = lammps_tools::data_field_int;
	if( d.df->type() != lammps_tools::data_field::INT ){
		return nullptr;
	}else{
		const dfi* dd = static_cast<const dfi*>( d.df );
		return dd->get_data_ptr();
	}
}

int lt_data_type( lt_data_field_handle d )
{
	return d.df->type();
}

const char *lt_data_name( lt_data_field_handle d )
{
	return d.df->name.c_str();
}


const std::vector<double> &lt_data_as_double_vec ( lt_data_field_handle d )
{
	return lammps_tools::data_as<double>( d.df );
}

const std::vector<int> &lt_data_as_int_vec ( lt_data_field_handle d )
{
	return lammps_tools::data_as<int>( d.df );
}
