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
enum LT_DATA_FIELD_TYPES {
	DATA_FIELD_DOUBLE = lammps_tools::data_field::DOUBLE,
	DATA_FIELD_INT	  = lammps_tools::data_field::INT
};


struct lt_data_field_handle
{
private:
	lammps_tools::data_field *df_;
	const lammps_tools::data_field *df;
	lammps_tools::data_field *df_rw;
	bool clean_up;
public:
	lt_data_field_handle() :
		df_(nullptr), df(df_), df_rw(nullptr), clean_up(false) {}

	explicit lt_data_field_handle( const lammps_tools::data_field *df );
	explicit lt_data_field_handle( lammps_tools::data_field *df );

	lt_data_field_handle( const char *name, int type, int size );
	~lt_data_field_handle();

	const lammps_tools::data_field *get() const { return df; }
	lammps_tools::data_field *get() { return df_rw; }
	bool is_rw() const
	{
		return df_rw != nullptr;
	}

	void set( const lammps_tools::data_field *field )
	{
		df = field;
		df_rw = nullptr;
	};
	void set_rw( lammps_tools::data_field *field )
	{
		df = field;
		df_rw = field;
	};
};

int           lt_data_field_size( const lt_data_field_handle *d );
const double *lt_data_as_double ( const lt_data_field_handle *d );
const int    *lt_data_as_int    ( const lt_data_field_handle *d );
int           lt_data_field_type( const lt_data_field_handle *d );
const char   *lt_data_field_name( const lt_data_field_handle *d );

double *lt_data_as_double_rw( lt_data_field_handle *d );
int    *lt_data_as_int_rw   ( lt_data_field_handle *d );


void lt_data_field_set_size( lt_data_field_handle *d, int size );
void lt_data_field_set_name( lt_data_field_handle *d, const char *name );


int lt_data_field_get_indexed_int_data( const lt_data_field_handle *d,
                                        int index );
int lt_data_field_set_indexed_int_data( lt_data_field_handle *d,
                                        int index, int val );

double lt_data_field_get_indexed_double_data( const lt_data_field_handle *d,
                                              int index );
int lt_data_field_set_indexed_double_data( lt_data_field_handle *d,
                                           int index, double val );



lt_data_field_handle *lt_new_data_field( const char *name, int dtype, int size );
lt_data_field_handle *lt_new_empty_data_field();

void lt_delete_data_field( lt_data_field_handle *d );

} // extern "C"

const std::vector<double> &lt_data_as_double_vec ( const lt_data_field_handle *d );
const std::vector<int>    &lt_data_as_int_vec ( const lt_data_field_handle *d );


// Ease of use template:
template <typename T> inline
const T *lt_data_as( const lt_data_field_handle *d )
{
	if     ( typeid(T) == typeid(double) ) return lt_data_as_double( d );
	else if( typeid(T) == typeid(double) ) return lt_data_as_int( d );

	return nullptr;
}



#endif // LT_DATA_FIELD_H
