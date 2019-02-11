#include "atom_type_info.hpp"
#include "block_data.hpp"
#include "domain.hpp"
#include "id_map.hpp"
#include "my_assert.hpp"

#include <iostream>

using namespace lammps_tools;

namespace lammps_tools {

block_data::block_data()
	: tstep( 0 ), N( 0 ), N_ghost(0), N_true(0),
	  N_types( 1 ), atom_style( ATOM_STYLE_ATOMIC ),
	  dom(), top(), ati(N_types), data(),
	  special_fields_by_name ( N_SPECIAL_FIELDS, "" ),
	  special_fields_by_index( N_SPECIAL_FIELDS, -1 ),
	  field_to_special_field_type(0)
{ }

block_data::block_data( std::size_t n_atoms )
	: tstep( 0 ), N( n_atoms ), N_ghost(0), N_true(0),
	  N_types( 1 ), atom_style( ATOM_STYLE_ATOMIC ),
	  dom(), top(), ati(N_types), data(),
	  special_fields_by_name ( N_SPECIAL_FIELDS, "" ),
	  special_fields_by_index( N_SPECIAL_FIELDS, -1 ),
	  field_to_special_field_type(0)
{ }

block_data::block_data( const block_data &o )
	: tstep( o.tstep ),
	  N( o.N ),
	  N_ghost( o.N_ghost ),
	  N_true( o.N_true ),
	  N_types( o.N_types ),
	  atom_style( o.atom_style ),
	  dom( o.dom ),
	  top( o.top ),
	  ati( o.ati ),
	  data( o.n_data_fields() ),
	  special_fields_by_name ( N_SPECIAL_FIELDS, "" ),
	  special_fields_by_index( N_SPECIAL_FIELDS, -1 ),
	  field_to_special_field_type( o.n_data_fields() )
{
	my_assert( __FILE__, __LINE__, data.size() == o.n_data_fields(),
	           "Data size mismatch after copy!" );

	for( std::size_t i = 0; i < o.n_data_fields(); ++i ){
		data[i] = copy( o.get_data( o[i].name) );
		field_to_special_field_type[i] = UNKNOWN;
	}

	for( int i = block_data::ID; i < block_data::N_SPECIAL_FIELDS; ++i ){
		const data_field *df = o.get_special_field(i);
		if( !df ) continue;

		special_fields_by_name[i] = df->name;

		for( std::size_t idx = 0; idx < data.size(); ++idx ){
			const data_field *df2 = data[idx];
			if( df2->name == special_fields_by_name[i] ){
				special_fields_by_index[i] = idx;
				field_to_special_field_type[idx] = i;
			}
		}
	}
}

block_data::~block_data()
{
	for( data_field *df : data ){
		delete df;
	}
	/*
	std::cerr << "Deleted block_data at " << this
	          << ", hope you don't need it anymore...\n";
	*/
}

data_field *block_data::get_data_rw( const std::string &name )
{
	for( data_field *df : data ){
		if( df->name == name ){
			return df;
		}
	}
	return nullptr;
}

const data_field *block_data::get_data( const std::string &name ) const
{
	for( const data_field *df : data ){
		if( df->name == name ){
			return df;
		}
	}

	return nullptr;
}


void block_data::add_field( const data_field &data_f, int special_field )
{
	my_assert( __FILE__, __LINE__,
	           data_f.size() == static_cast<std::size_t>(N),
	           "Atom number mismatch on add_field! Call set_natoms first!");
	// Check if this name is already in block or not.
	if( get_data( data_f.name ) != nullptr ){
		my_runtime_error( __FILE__, __LINE__,
		                  "Named data already in block_data" );
	}
	// Now you need to copy the data.
	data_field *cp = copy( &data_f );
	int index = data.size();
	data.push_back( cp );

	field_to_special_field_type.push_back( UNKNOWN );

	// std::cerr << "Now block_data has " << data.size() << " data fields.\n";
	// Ignore some keys that are not unique for example:
	if( !is_legal_special_field( special_field ) ){
		// print_internal_state();
		return;
	}
	my_assert( __FILE__, __LINE__,
	           special_fields_by_index[special_field] == -1,
	           "Special field already set!" );

	special_fields_by_name[special_field]  = data_f.name;
	special_fields_by_index[special_field] = index;
	field_to_special_field_type[index] = special_field;

	// If the field is molecule IDs, switch atom style to molecular:
	if( special_field == block_data::special_fields::MOL ){
		atom_style = ATOM_STYLE_MOLECULAR;
	}

	//print_internal_state();
}

void block_data::print_internal_state()
{
	std::cerr << " *** my " << data.size() << " fields are:";
	for( const data_field *df : data ){
		std::cerr << " " << df->name;
	}
	std::cerr << "\n";
	std::cerr << " *** special_fields_by_name is\n        ";
	for( const std::string &s : special_fields_by_name ){
		std::cerr << " " << s;
	}
	std::cerr << "\n *** special_fields_by_index is\n        ";
	for( int i : special_fields_by_index ){
		std::cerr << " " << i;
	}
	std::cerr << "\n *** field_to_special_field_type is\n     ";
	for( std::size_t i = 0; i < data.size(); ++i ){
		std::cerr << " " << field_to_special_field_type[i];
	}
	std::cerr << "\n\n";
}

block_data &block_data::operator=( block_data o )
{
	if( this != &o ){
		using std::swap;
		block_data n(o);
		swap( *this, n );
	}
	return *this;
}



std::size_t block_data::n_data_fields() const
{
	return data.size();
}


void block_data::set_natoms( std::size_t new_size )
{
	for( data_field *d : data ){
		my_assert( __FILE__, __LINE__, d,
		           "Data field not initialised in attempted resize!" );
		d->resize(new_size);
	}
	N = new_size;
}


void block_data::set_ntypes( std::size_t N )
{
	N_types = N;
	ati.set_size( N_types );
}


void block_data::set_special_field( const std::string &name, int field )
{
	std::size_t index = 0;
	my_assert( __FILE__, __LINE__, is_legal_special_field( field ),
	           "Invalid field in set_special_field!" );
	while( index < data.size() && data[index]->name != name ){
		++index;
	}
	if( index < data.size() ){
		special_fields_by_name[field]  = name;
		special_fields_by_index[field] = index;
	}
}

int block_data::get_special_field_type( const std::string& name ) const
{
	std::size_t idx = name2index( name );
	return get_special_field_type( idx );
}

std::size_t block_data::name2index( const std::string &name ) const
{
	std::size_t index = 0;
	while( index < data.size() && data[index]->name != name ) ++index;
	return index;
}

int block_data::get_special_field_type( int idx ) const
{
	if( idx >= field_to_special_field_type.size() ){
		return UNKNOWN;
	}
	return field_to_special_field_type[ idx ];
}




std::string block_data::get_special_field_name( int field ) const
{
	my_assert( __FILE__, __LINE__, is_legal_special_field( field ),
	           "Invalid field in get_special_field_name!" );

	return special_fields_by_name[field];
}




data_field *block_data::remove_field( const std::string &name,
                                      int &special_field )
{
	// Delete the ptr from all vectors and maps.
	int index = 0;
	data_field *df = nullptr;
	for( data_field *d : data ){
		if( d->name == name ){
			df = d;
			break;
		}
		++index;
	}
	if( !df ) return nullptr;
	data.erase( data.begin() + index );
	//field_to_special_field_type.erase(
	//	field_to_special_field_type.begin() + index );

	bool found = false;
	for( std::size_t i = 0; i < special_fields_by_name.size(); ++i ){
		if( special_fields_by_name[i] == name ){
			special_field = i;
			found = true;
			break;
		}
	}
	if( found ){
		special_fields_by_name[ special_field ] = "";
		special_fields_by_index[ special_field ] = -1;

		for( int &i : special_fields_by_index ){
			if( i > index ){
				--i;
			}
		}
	}

	// Check if this leaves the internal state correct:
	//print_internal_state();
	return df;
}


/// Get read/write pointer to special data field of given kind.
data_field *block_data::get_special_field_rw( int field )
{
	if( special_fields_by_index[field] != -1 ){
		return data[ special_fields_by_index[field] ];
	}else{
		return nullptr;
	}
}

/// Get read-only pointer to special data field of given kind.
const data_field *block_data::get_special_field( int field ) const
{
	if( special_fields_by_index[field] != -1 ){
		int entry = special_fields_by_index[field];
		return data[ entry ];
	}else{
		return nullptr;
	}
}



int block_data::n_special_fields() const
{
	int c = 0;
	for( int i : special_fields_by_index ){
		if( i >= 0 ) ++c;
	}
	return c;
}


const data_field &block_data::operator[]( int i ) const
{
	return *data[i];
}

data_field &block_data::operator[]( int i )
{
	return *data[i];
}


data_field *block_data::get_data_rw( int index )
{
	return data[index];
}


// ******************   Non-member functions:    ************************
void swap( block_data &f, block_data &s )
{
	// using std::swap;

	f.tstep = s.tstep;
	f.N = s.N;
	f.N_ghost = s.N_ghost;
	f.N_true = s.N_true;
	f.N_types = s.N_types;
	f.atom_style = s.atom_style;
	swap( f.dom, s.dom );
	swap( f.top, s.top );
	swap( f.data, s.data );
	swap( f.special_fields_by_name, s.special_fields_by_name );
	swap( f.special_fields_by_index, s.special_fields_by_index );
	swap( f.ati, s.ati );
	swap( f.field_to_special_field_type, s.field_to_special_field_type );
}



bool grab_common_fields( const block_data &b,
                         const std::vector<std::string> &fields,
                         std::vector<int> &id, std::vector<int> &type,
                         std::vector<double> &x, std::vector<double> &y,
                         std::vector<double> &z )
{
	using dfd = data_field_double;
	using dfi = data_field_int;

	my_assert( __FILE__, __LINE__, fields.size() == 5,
	           "fields has incorrect size!" );

	const data_field *df_i = b.get_data( fields[0] );
	const data_field *df_t = b.get_data( fields[1] );
	const data_field *df_x = b.get_data( fields[2] );
	const data_field *df_y = b.get_data( fields[3] );
	const data_field *df_z = b.get_data( fields[4] );

	if( !(df_i && df_t && df_x && df_y && df_z) ){
		return false;
	}

	id   = static_cast<const dfi*>(df_i)->get_data();
	type = static_cast<const dfi*>(df_t)->get_data();
	x    = static_cast<const dfd*>(df_x)->get_data();
	y    = static_cast<const dfd*>(df_y)->get_data();
	z    = static_cast<const dfd*>(df_z)->get_data();

	return true;
}




bool is_special_field_int( int special_field )
{
	switch( special_field ){
		default:
		case block_data::UNKNOWN:
		case block_data::X:
		case block_data::Y:
		case block_data::Z:
		case block_data::VX:
		case block_data::VY:
		case block_data::VZ:
			return false;

		case block_data::ID:
		case block_data::TYPE:
		case block_data::MOL:
		case block_data::IX:
		case block_data::IY:
		case block_data::IZ:
			return true;
	}

}



bool is_legal_special_field( int special_field )
{
	return (special_field >= block_data::ID) &&
		(special_field < block_data::N_SPECIAL_FIELDS);
}



block_data filter_block( const block_data &b, const std::vector<int> &ids )
{
	block_data new_block( b );
	bigint new_size = ids.size();
	new_block.set_natoms( new_size );
	std::size_t nfields = b.n_data_fields();

	// Loop over all data fields and extract new ones from them,
	// then replace the one in new_block with those.
	int field_id = block_data::special_fields::ID;
	const data_field_int *df_ids = static_cast<const data_field_int*>(
		b.get_special_field( field_id ) );

	id_map im( df_ids->get_data() );
	new_block.clear();

	for( std::size_t i = 0; i < nfields; ++i ){
		const data_field *df_i = &b[i];
		const std::string &field_name = df_i->name;
		int type = df_i->type();

		data_field_int *dfi;
		data_field_double *dfd;

		const data_field_int *dfi_old = nullptr;
		const data_field_double *dfd_old = nullptr;

		data_field *df_new;

		if( type == data_field::INT ){
			dfi = new data_field_int( field_name, new_size );
			dfi_old = static_cast<const data_field_int*>( df_i );
			dfd = nullptr;
			df_new = dfi;
		}else{
			dfd = new data_field_double( field_name, new_size );
			dfd_old = static_cast<const data_field_double*>( df_i );
			dfi = nullptr;
			df_new = dfd;
		}
		df_new->name = df_i->name;
		int special_field_type = b.get_special_field_type( i );

		std::size_t k = 0;
		for( std::size_t j : ids ){
			std::size_t idx = im[j];
			if( dfi ){
				(*dfi)[k] = (*dfi_old)[idx];
				++k;
				continue;
			}else if( dfd ){
				(*dfd)[k] = (*dfd_old)[idx];
				++k;
				continue;
			}

			my_runtime_error( __FILE__, __LINE__,
			                  "data_field ptr assignment failed!" );
		}
		if( dfi ){
			new_block.add_field( *dfi, special_field_type );
			delete dfi;
		}else if( dfd ){
			new_block.add_field( *dfd, special_field_type );
			delete dfd;
		}
	}



	return new_block;
}



void block_data::clear()
{
	while( n_data_fields() > 0 ){
		int special_field_type = 0;

		data_field *df = remove_field( data[0]->name,
		                               special_field_type );
		delete df;
	}
}


void block_data::add_ghost_atoms( double rc, int dims )
{
	// For each atom, if it is within rc of the edge, create a copy
	// for each periodic boundary it is within rc of.
	my_assert( __FILE__, __LINE__, !have_ghost_atoms(),
	           "Adding ghost atoms more than once not supported!" );
	std::cerr << "Adding ghost atoms...\n";
	int ghost_added = 0;

	using dfd = data_field_double;
	dfd *xx = static_cast<dfd*>(get_special_field_rw( block_data::X ));
	dfd *yy = static_cast<dfd*>(get_special_field_rw( block_data::Y ));
	dfd *zz = static_cast<dfd*>(get_special_field_rw( block_data::Z ));

	bool periodic_x = dom.periodic & domain::BIT_X;
	bool periodic_y = dom.periodic & domain::BIT_Y;
	bool periodic_z = dom.periodic & domain::BIT_Z;
	if (dims == 2) periodic_z = false;

	bigint old_N = N;

	for( bigint i = 0; i < old_N; ++i ){
		double x = (*xx)[i];
		double y = (*yy)[i];
		double z = (*zz)[i];

		int flip_x = 0; // 1 means "add to the right" with right being
		int flip_y = 0; // at xhi, -1 means "add to the left" with left
		int flip_z = 0; // being at xlo, 0 means nothing needed.

		if (periodic_x) {
			// Suppose xlo = 0, xhi = 10, rc = 1.
			// Then any position [0,1] and [9,10] needs
			// mirroring. That means [.2 - 0 ] < rc
			// or [10 - 9.8 ] < rc:
			if     (x - dom.xlo[0] <= rc) flip_x =  1;
			else if(dom.xhi[0] - x <  rc) flip_x = -1;
		}
		if (periodic_y) {
			if     (y - dom.xlo[1] <= rc) flip_y =  1;
			else if(dom.xhi[1] - y <  rc) flip_y = -1;
		}
		if (periodic_z) {
			if     (z - dom.xlo[2] <= rc) flip_z =  1;
			else if(dom.xhi[2] - z <  rc) flip_z = -1;
		}

		if (!flip_x && !flip_y && !flip_z) continue;

		// Handle all the cases:
		double Lx = dom.xhi[0] - dom.xlo[0];
		double Ly = dom.xhi[1] - dom.xlo[1];
		double Lz = dom.xhi[2] - dom.xlo[2];

		double dx = x + flip_x*Lx;
		double dy = y + flip_y*Ly;
		double dz = z + flip_z*Lz;


		if (flip_x) {
			bigint new_idx = clone_particle(i);
			(*xx)[new_idx] = dx;
			++ghost_added;
		}
		if (flip_y) {
			bigint new_idx = clone_particle(i);
			(*yy)[new_idx] = dy;
			++ghost_added;
		}
		if (flip_z) {
			bigint new_idx = clone_particle(i);
			(*zz)[new_idx] = dz;
			++ghost_added;
		}

		if (flip_x && flip_y) {
			bigint new_idx = clone_particle(i);

			(*xx)[new_idx] = dx;
			(*yy)[new_idx] = dy;
			++ghost_added;
		}

		if (flip_x && flip_z) {
			bigint new_idx = clone_particle(i);
			(*xx)[new_idx] = dx;
			(*zz)[new_idx] = dz;
			++ghost_added;
		}

		if (flip_y && flip_z) {
			bigint new_idx = clone_particle(i);
			(*yy)[new_idx] = dy;
			(*zz)[new_idx] = dz;
			++ghost_added;
		}

		if (flip_x && flip_y && flip_z) {
			bigint new_idx = clone_particle(i);
			(*xx)[new_idx] = dx;
			(*yy)[new_idx] = dy;
			(*zz)[new_idx] = dz;
			++ghost_added;
		}

	}

	std::cerr << "Added " << ghost_added << " ghost atoms.\n";
	bigint N_added = N - old_N;
	my_assert(__FILE__, __LINE__, N_added == ghost_added,
	          "Mismatch in number of ghosts added!");
	N_ghost = ghost_added;
	N_true = old_N;
	my_assert(__FILE__, __LINE__, N = N_true + N_ghost,
	          "Mismatch in number of ghosts/true!");

	// Extend the domain:
	dom.xlo[0] -= rc;
	dom.xhi[0] += rc;
	dom.xlo[1] -= rc;
	dom.xhi[1] += rc;
	if (dims == 3) {
		dom.xlo[2] -= rc;
		dom.xhi[2] += rc;
	}
}



bigint block_data::clone_particle( bigint idx )
{
	my_assert(__FILE__, __LINE__, idx >= 0 && idx < N, "Index out of range!");
	std::size_t new_pos = N;
	set_natoms(N+1);
	using dfd = data_field_double;
	using dfi = data_field_int;

	int n_props = 0;

	for ( data_field *d : data ) {
		if (d->type() == data_field::DOUBLE) {
			double val = (*static_cast<dfd*>(d))[idx];
			(*static_cast<dfd*>(d))[new_pos] = val;
		}else {
			int val = (*static_cast<dfi*>(d))[idx];
			(*static_cast<dfi*>(d))[new_pos] = val;
		}
		++n_props;
	}
	return new_pos;
}



void block_data::unwrap_image(std::size_t idx, double dest[3]) const
{
	const std::vector<double> &x = data_as<double>(get_special_field(block_data::X));
	const std::vector<double> &y = data_as<double>(get_special_field(block_data::Y));
	const std::vector<double> &z = data_as<double>(get_special_field(block_data::Z));

	const std::vector<int> &ix = data_as<int>(get_special_field(block_data::IX));
	const std::vector<int> &iy = data_as<int>(get_special_field(block_data::IY));
	const std::vector<int> &iz = data_as<int>(get_special_field(block_data::IZ));

	dest[0] = x[idx];
	dest[1] = y[idx];
	dest[2] = z[idx];
	int flags[3] = {ix[idx], iy[idx], iz[idx]};
	dom.unwrap_image(dest, flags);
}



} // namespace lammps_tools
