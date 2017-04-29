#include "data_field.hpp"

using namespace lammps_tools;

namespace lammps_tools {

void swap( data_field &f, data_field &s )
{
	// This is painful... Not really.
	int type_f = f.type();
	int type_s = s.type();
	my_assert( __FILE__, __LINE__, type_f == type_s,
	           "Cannot swap different types!" );

	switch( type_f ){
		default: {
			my_runtime_error( __FILE__, __LINE__, "Unkown type in swap!" );
			break;
		}
		case data_field::INT:{
			data_field_int &d_f = dynamic_cast<data_field_int&>(f);
			data_field_int &d_s = dynamic_cast<data_field_int&>(s);
			swap( d_f, d_s );
		}
		case data_field::DOUBLE:{
			data_field_double &d_f = dynamic_cast<data_field_double&>(f);
			data_field_double &d_s = dynamic_cast<data_field_double&>(s);
			swap( d_f, d_s );			
		}
	}
}

} // namespace lammps_tools
