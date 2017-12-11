#include "lt_block_data.h"
#include "lt_data_field.h"

// #include "../cpp_lib/block_data_access.hpp"

#include "../cpp_lib/data_field.hpp"


lt_block_data_handle *lt_new_block_data_handle()
{
	lt_block_data_handle *bdh = new lt_block_data_handle;
	return bdh;
}

void lt_delete_block_data_handle( lt_block_data_handle *bdh )
{
	delete bdh;
}

bool lt_has_special_field( lt_block_data_handle *bdh, int special_field )
{
	// Cannot throw:
	return bdh->bd->get_special_field( special_field ) != nullptr;
}


const std::vector<double> &lt_special_field_double( lt_block_data_handle *bdh,
                                                    int special_field )
{
	lammps_tools::data_field *df
		= bdh->bd->get_special_field_rw( special_field );
	if( !df ){
		std::cerr << "No data field for special field type "
		          << special_field << "!\n";
		std::terminate();
	}

	lt_data_field_handle dfh;
	dfh.set(df);
	return lt_data_as_double_vec( &dfh );
}


const std::vector<int> &lt_special_field_int( lt_block_data_handle *bdh,
                                              int special_field )
{
	lammps_tools::data_field *df
		= bdh->bd->get_special_field_rw( special_field );
	if( !df ){
		std::cerr << "No data field for special field type "
		          << special_field << "!\n";
		std::terminate();
	}
	lt_data_field_handle dfh;
	dfh.set( df );
	return lt_data_as_int_vec( &dfh );
}


lt_data_field_handle lt_data_by_name( lt_block_data_handle *bdh, const char *name )
{
	const std::string n(name);
	lt_data_field_handle ldf;
	ldf.set(nullptr);
	const lammps_tools::block_data &b = *bdh->bd;

	for( std::size_t i = 0; i < b.n_data_fields(); ++i ){
		if( n == b[i].name ){
			const lammps_tools::data_field *df = &b[i];
			ldf.set( df );
			break;
		}
	}
	return ldf;
}


lt_data_field_handle lt_data_by_index( lt_block_data_handle *bdh, int i )
{
	lt_data_field_handle ldf;
	ldf.set(nullptr);
	const lammps_tools::block_data &b = *bdh->bd;
	std::size_t ii = i;
	if( i < 0 || ii >= b.n_data_fields() ){
		std::cerr << "Index " << i << " is out of range! "
		          << "Ignoring call to lt_data_by_index!\n";
		return ldf;
	}

	// This cannot go wrong as the index is already asserted to be OK.
	ldf.set(&b[i]);
	return ldf;
}


int lt_n_data_fields( lt_block_data_handle *bdh )
{
	return bdh->bd->n_data_fields();
}

void lt_block_data_add_data_field( lt_block_data_handle *bdh,
                                   const lt_data_field_handle *dfh )
{
	bdh->bd->add_field( *dfh->get() );
}

void lt_block_data_add_special_field( lt_block_data_handle *bdh,
                                      const lt_data_field_handle *dfh,
                                      int type )
{
	bdh->bd->add_field( *dfh->get(), type );
}


void lt_block_data_set_meta( lt_block_data_handle *bdh,
                             lammps_tools::bigint tstep,
                             lammps_tools::bigint natoms,
                             double xlo, double xhi,
                             double ylo, double yhi,
                             double zlo, double zhi,
                             int periodic_bits, int atom_style )
{
	bdh->bd->tstep  = tstep;
	bdh->bd->set_natoms( natoms );
	bdh->bd->atom_style = atom_style;

	lt_block_data_set_domain( bdh, xlo, xhi, ylo, yhi, zlo, zhi,periodic_bits );
}

void lt_block_data_set_domain( lt_block_data_handle *bdh,
                               double xlo, double xhi,
                               double ylo, double yhi,
                               double zlo, double zhi,
                               int periodic_bits )
{
	bdh->bd->dom.xlo[0] = xlo;
	bdh->bd->dom.xlo[1] = ylo;
	bdh->bd->dom.xlo[2] = zlo;

	bdh->bd->dom.xhi[0] = xhi;
	bdh->bd->dom.xhi[1] = yhi;
	bdh->bd->dom.xhi[2] = zhi;

	bdh->bd->dom.periodic = periodic_bits;
}


template <typename T>
void lt_block_data_set_data_impl( lt_block_data_handle *bdh,
                                  const char *name,
                                  const std::vector<T> &data )
{
	lt_data_field_handle dfh = lt_data_by_name( bdh, name );
	if( dfh.get() == nullptr ){
		std::cerr << "Could not find field named " << name << "!\n";
		return;
	}

	typedef typename std::conditional<std::is_same<T, double>::value,
	                                  lammps_tools::data_field_double,
	                                  lammps_tools::data_field_int>::type df_type;
	std::cerr << "Trying to find " << name << " among fields.\n";
	std::cerr << "Trying to cast " << dfh.get() << " to "
	          << typeid(df_type).name()
	          << " through " << dfh.get() << "\n";
	df_type *field = static_cast<df_type *>( dfh.get() );

	if( data.size() != field->size() ){
		std::cerr << "Data size mismatch for field " << name << "!\n";
		return;
	}

	std::cerr << "Field is " << field << ", dfh = " << dfh.get() << "\n";
	for( std::size_t i = 0; i < data.size(); ++i ){
		(*field)[i] = data[i];
	}

}

/**
   \brief Overwrites named data assuming it's doubles.
*/
void lt_block_data_set_data_double( lt_block_data_handle *bdh,
                                    const char *name,
                                    const std::vector<double> &data )
{
	lt_block_data_set_data_impl<double>( bdh, name, data );
}

/**
   \brief Overwrites named data assuming it's ints.
*/
void lt_block_data_set_data_int( lt_block_data_handle *bdh,
                                 const char *name,
                                 const std::vector<int> &data )
{
	lt_block_data_set_data_impl<int>( bdh, name, data );
}


/**
   \brief Prints some stats about the block data.
*/
void lt_block_data_print_stats( lt_block_data_handle *bdh )
{
	lammps_tools::block_data *b = bdh->bd;
	std::cerr << "Block data handle at " << bdh << ":\n";
	std::cerr << "  Block data at " << b << "\n";
	std::cerr << "  Atom stype = " << b->atom_style << "\n";
	std::cerr << "  t = " << bdh->time_step() << "\n";
	std::cerr << "  " << bdh->n_atoms() << " atoms\n";
	std::cerr << "  " << bdh->n_types() << " types\n";
	std::cerr << "  " << b->n_data_fields() << " data fields:\n";

	for( uint i = 0; i < b->n_data_fields(); ++i ){
		std::cerr << "    " << (*b)[i].name << "\n";
	}
	std::cerr << "  " << b->n_special_fields() << " special fields:\n";

	for( int i = 0; i < lammps_tools::block_data::N_SPECIAL_FIELDS; ++i ){
		if( lt_has_special_field( bdh, i ) ){
			const lammps_tools::data_field *sf = b->get_special_field(i);
			std::cerr << "    " << sf->name << "\n";
		}
	}


	std::cerr << "\n";
}


void lt_block_data_remove_field( lt_block_data_handle *bdh, const char *name )
{
	int spec_type = lammps_tools::block_data::UNKNOWN;
	lammps_tools::data_field *df = bdh->bd->remove_field( name, spec_type );
	if( !df ){
		std::cerr << "Warning: Data field named " << name
		          << " could not be found!\n";
		return;
	}
	delete df;
}

void lt_block_data_swap_fields( lt_block_data_handle *bdh, const char *name,
                                const lt_data_field_handle *new_df )
{
	int spec_type = lammps_tools::block_data::special_fields::UNKNOWN;
	lammps_tools::data_field *df = bdh->bd->remove_field( name, spec_type );

	if( spec_type != lammps_tools::block_data::special_fields::UNKNOWN ){
		lt_block_data_add_special_field( bdh, new_df, spec_type );
	}else{
		lt_block_data_add_data_field( bdh, new_df );
	}
	delete df;
}

const std::vector<double>
lt_block_data_get_domain_xlo_vec( const lt_block_data_handle *bdh )
{
	std::vector<double> xlo(3);
	xlo[0] = bdh->bd->dom.xlo[0];
	xlo[1] = bdh->bd->dom.xlo[1];
	xlo[2] = bdh->bd->dom.xlo[2];
	return xlo;
}

const std::vector<double>
lt_block_data_get_domain_xhi_vec( const lt_block_data_handle *bdh )
{
	std::vector<double> xhi(3);
	xhi[0] = bdh->bd->dom.xhi[0];
	xhi[1] = bdh->bd->dom.xhi[1];
	xhi[2] = bdh->bd->dom.xhi[2];
	return xhi;
}

const double *lt_block_data_get_domain_xlo( const lt_block_data_handle *bdh )
{
	return bdh->bd->dom.xlo;
}

const double *lt_block_data_get_domain_xhi( const lt_block_data_handle *bdh )
{
	return bdh->bd->dom.xhi;
}

int lt_block_data_get_domain_periodic( const lt_block_data_handle *bdh )
{
	return bdh->bd->dom.periodic;
}


void lt_block_data_filter( lt_block_data_handle *dest, int size, const void *ids,
                           const lt_block_data_handle *src )
{
	const int *id = static_cast<const int*>( ids );
	std::vector<int> id_vec( id, id + size );
	// Make sure the old storage is properly deleted.

	if( dest->bd ){
		delete dest->bd;
	}
	lammps_tools::block_data temp_b =
		lammps_tools::filter_block( src->get_const_ref(), id_vec );
	dest->bd = new lammps_tools::block_data( temp_b );
}
