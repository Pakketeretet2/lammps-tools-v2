#include "lt_data_field.h"

#include "../cpp_lib/data_field.hpp"

lt_data_field_handle::lt_data_field_handle( const char *name,
                                            int type, int size )
	: df_(nullptr), df(df_), df_rw(nullptr), clean_up(false)
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
	df_rw = df_;
	clean_up = true;
}

lt_data_field_handle::lt_data_field_handle( const lammps_tools::data_field *df )
	: df_(nullptr), df(df), df_rw(nullptr), clean_up(false)
{}

lt_data_field_handle::lt_data_field_handle( lammps_tools::data_field *df )
	: df_(nullptr), df(df), df_rw(df), clean_up(false)
{}

lt_data_field_handle::~lt_data_field_handle()
{
	// Leave df dangling, but delete df_ if it was allocated:
	if( clean_up && df_ ){
		/*
		std::cerr << "This is lt_data_field_handle at " << this
		          << ", cleaning up data_field at " << df_ << "\n";
		*/
		delete df_;
	}
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
	return d->get()->size();
}

const double *lt_data_as_double( const lt_data_field_handle *d )
{
	using dfd = lammps_tools::data_field_double;
	if( d->get()->type() != lammps_tools::data_field::DOUBLE ){
		return nullptr;
	}else{
		const dfd* dd = static_cast<const dfd*>( d->get() );
		return dd->get_data_ptr();
	}
}

const int *lt_data_as_int( const lt_data_field_handle *d )
{
	using dfi = lammps_tools::data_field_int;
	if( d->get()->type() != lammps_tools::data_field::INT ){
		return nullptr;
	}else{
		const dfi* dd = static_cast<const dfi*>( d->get() );
		return dd->get_data_ptr();
	}
}

int lt_data_field_type( const lt_data_field_handle *d )
{
	return d->get()->type();
}

const char *lt_data_field_name( const lt_data_field_handle *d )
{
	return d->get()->name.c_str();
}


const std::vector<double> &lt_data_as_double_vec ( const lt_data_field_handle *d )
{
	try{
		const std::vector<double> &dvec =
			lammps_tools::data_as<double>( d->get() );
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
			lammps_tools::data_as<int>( d->get() );
		return dvec;
	}catch( std::runtime_error &e ){
		std::cerr << "Error occured in lt_data_as_int_vec! "
		          << e.what() << "\n";
		std::terminate();
	}
}


void lt_data_field_set_size( lt_data_field_handle *d, int size )
{
	d->get()->resize( size );
}

void lt_data_field_set_name( lt_data_field_handle *d, const char *name )
{
	d->get()->name = name;
}



template <typename T>
int lt_data_field_set_indexed_data( lt_data_field_handle *d,
                                    int index, T val )
{
	if( !d->is_rw() ){
		std::cerr << "Cannot mutate data in const data_field!\n";
		return -1;
	}

	lammps_tools::data_field *df = d->get();
	if( !df ){
		std::cerr << "Failed to grab data from handle @ " << d << "!\n";
		return -2;
	}

	constexpr const int TYPE = lammps_tools::data_field_type_id<T>::TYPE;
	typedef lammps_tools::data_field_der<T, TYPE> der_type;
	der_type *dfd = static_cast<der_type*> ( df );
	dfd->set_val_by_index( index, val );
	return 0;
}

template <typename T>
T lt_data_field_get_indexed_data( const lt_data_field_handle *d, int index )
{
	return lammps_tools::data_as<T>( d->get() )[index];
}

int lt_data_field_get_indexed_int_data( const lt_data_field_handle *d,
                                        int index )
{
	return lt_data_field_get_indexed_data<int>( d, index );
}


int lt_data_field_set_indexed_int_data( lt_data_field_handle *d,
                                        int index, int val )
{
	return lt_data_field_set_indexed_data<int>( d, index, val );
}

double lt_data_field_get_indexed_double_data( const lt_data_field_handle *d,
                                              int index )
{
	return lt_data_field_get_indexed_data<double>( d, index );
}

int lt_data_field_set_indexed_double_data( lt_data_field_handle *d,
                                            int index, double val )
{
	return lt_data_field_set_indexed_data<double>( d, index, val );
}

double *lt_data_as_double_rw( lt_data_field_handle *d )
{
	using dfd = lammps_tools::data_field_double;
	if( d->get()->type() != lammps_tools::data_field::DOUBLE ){
		return nullptr;
	}else{
		dfd* dd = static_cast<dfd*>( d->get() );
		std::vector<double> &vec = dd->get_data_rw();
		return vec.data();
	}
}

int *lt_data_as_int_rw( lt_data_field_handle *d )
{
	using dfi = lammps_tools::data_field_int;
	if( d->get()->type() != lammps_tools::data_field::INT ){
		return nullptr;
	}else{
		dfi* dd = static_cast<dfi*>( d->get() );
		std::vector<int> &vec = dd->get_data_rw();
		return vec.data();
	}
}
