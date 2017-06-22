#include "lt_data_field.h"

#include "../cpp_lib/data_field.hpp"

lt_data_field_handle::lt_data_field_handle( const char *name,
                                            int type, int size )
	: df(df_), df_(nullptr)
{
	switch( type ){
		case DATA_FIELD_DOUBLE:
			df_ = new lammps_tools::data_field_double( name, size );
			break;
		case DATA_FIELD_INT:
			df_ = new lammps_tools::data_field_int( name, size );
			break;
		default:
			// Error?
			std::cerr << "Unknown data type " << type
			          << "! Aborting!\n";
			std::terminate();
	}
	df = df_;
}

lt_data_field_handle::~lt_data_field_handle()
{
	// Leave df dangling, but delete df_ if it was allocated:
	if( df_ ) delete df_;
}


void lt_delete_data_field( lt_data_field_handle *d )
{
	delete d;
}


lt_data_field_handle *lt_new_data_field( const char *name, int dtype, int size )
{
	auto n = new lt_data_field_handle( name, dtype, size );
	return n;
}



int lt_data_field_size( const lt_data_field_handle *d )
{
	return d->df->size();
}

const double *lt_data_as_double( const lt_data_field_handle *d )
{
	using dfd = lammps_tools::data_field_double;
	if( d->df->type() != lammps_tools::data_field::DOUBLE ){
		return nullptr;
	}else{
		const dfd* dd = static_cast<const dfd*>( d->df );
		return dd->get_data_ptr();
	}
}

const int *lt_data_as_int( const lt_data_field_handle *d )
{
	using dfi = lammps_tools::data_field_int;
	if( d->df->type() != lammps_tools::data_field::INT ){
		return nullptr;
	}else{
		const dfi* dd = static_cast<const dfi*>( d->df );
		return dd->get_data_ptr();
	}
}

int lt_data_field_type( const lt_data_field_handle *d )
{
	return d->df->type();
}

const char *lt_data_field_name( const lt_data_field_handle *d )
{
	return d->df->name.c_str();
}


const std::vector<double> &lt_data_as_double_vec ( const lt_data_field_handle *d )
{
	try{
		const std::vector<double> &dvec =
			lammps_tools::data_as<double>( d->df );
		return dvec;
	}catch( std::runtime_error &e ){
		std::cerr << "Error occured in lt_data_as_double_vec! "
		          << e.what() << "!\n";
		std::terminate();
	}
}

const std::vector<int> &lt_data_as_int_vec ( const lt_data_field_handle *d )
{
	try{
		const std::vector<int> &dvec =
			lammps_tools::data_as<int>( d->df );
		return dvec;
	}catch( std::runtime_error &e ){
		std::cerr << "Error occured in lt_data_as_int_vec! "
		          << e.what() << "\n";
		std::terminate();
	}
}


void lt_data_field_set_size( lt_data_field_handle *d, int size )
{
	d->df_rw()->resize( size );
}

void lt_data_field_set_name( lt_data_field_handle *d, const char *name )
{
	d->df_rw()->name = name;
}


int lt_data_field_get_indexed_int_data( const lt_data_field_handle *d,
                                        int index )
{
	return lammps_tools::data_as<int>( d->df )[index];
}

void lt_data_field_set_indexed_int_data( lt_data_field_handle *d,
                                         int index, int val )
{
	lammps_tools::data_as_rw<int>( d->df_rw() )[index] = val;
}

double lt_data_field_get_indexed_double_data( const lt_data_field_handle *d,
                                              int index )
{
	return lammps_tools::data_as<double>( d->df )[index];
}

void lt_data_field_set_indexed_double_data( lt_data_field_handle *d,
                                            int index, double val )
{
	lammps_tools::data_as_rw<double>( d->df_rw() )[index] = val;
}
