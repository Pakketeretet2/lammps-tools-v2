#ifndef LT_DATA_FIELD_H
#define LT_DATA_FIELD_H

/**
   \file lt_data_field.h

   Exposes some features of data_fields through C interface.
*/

#include "../cpp_lib/data_field.hpp"

extern "C" {

/**
   \brief contains a pointer to an instantiated data field.
*/
struct lt_data_field_handle
{
	const lammps_tools::data_field *df;
};

int           lt_data_field_size( lt_data_field_handle d );
const double *lt_data_as_double ( lt_data_field_handle d );
const int    *lt_data_as_int    ( lt_data_field_handle d );
int           lt_data_type      ( lt_data_field_handle d );
const char   *lt_data_name      ( lt_data_field_handle d );


} // extern "C"

const std::vector<double> &lt_data_as_double_vec ( lt_data_field_handle d );
const std::vector<int> &lt_data_as_int_vec ( lt_data_field_handle d );


#endif // LT_DATA_FIELD_H
