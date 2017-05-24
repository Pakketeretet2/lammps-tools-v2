#include "block_data.hpp"
#include "my_assert.hpp"

using namespace lammps_tools;

namespace lammps_tools {

block_data::block_data()
	: tstep( 0 ), N( 0 ), N_types( 1 ), atom_style( ATOM_STYLE_ATOMIC ),
	  special_fields_by_name ( N_SPECIAL_FIELDS, "" ),
	  special_fields_by_index( N_SPECIAL_FIELDS, -1 )
{ }

block_data::block_data( std::size_t n_atoms )
	: tstep( 0 ), N( n_atoms ), N_types( 1 ), atom_style( ATOM_STYLE_ATOMIC ),
	  special_fields_by_name ( N_SPECIAL_FIELDS, "" ),
	  special_fields_by_index( N_SPECIAL_FIELDS, -1 )
{ }

block_data::block_data( const block_data &o )
	: tstep( o.tstep ),
	  N( o.N ),
	  N_types( o.N_types ),
	  atom_style( o.atom_style ),
	  dom( o.dom ),
	  top( o.top ),
	  data( o.n_data_fields() ),
	  special_fields_by_name ( N_SPECIAL_FIELDS, "" ),
	  special_fields_by_index( N_SPECIAL_FIELDS, -1 )
{
	my_assert( __FILE__, __LINE__, data.size() == o.n_data_fields(),
	           "Data size mismatch after copy!" );

	for( std::size_t i = 0; i < o.n_data_fields(); ++i ){
		data[i] = copy( o.get_data( o[i].name) );
	}

	for( int i = block_data::ID; i < block_data::N_SPECIAL_FIELDS; ++i ){
		const data_field *df = o.get_special_field(i);
		if( !df ) continue;

		special_fields_by_name[i] = df->name;

		for( int idx = 0; idx < data.size(); ++idx ){
			const data_field *df = data[idx];
			if( df->name == special_fields_by_name[i] ){
				special_fields_by_index[i] = idx;
			}
		}
	}
}

block_data::~block_data()
{
	for( data_field *df : data ){
		delete df;
	}
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


void block_data::add_field( const data_field &data_f, int special_field)
{
	my_assert( __FILE__, __LINE__, data_f.size() == N,
	           "Atom number mismatch on add_field! Call set_natoms first!");
	// Check if this name is already in block or not.
	if( get_data( data_f.name ) != nullptr ){
		std::cerr << "Name " << data_f.name << " was already "
		          << "in block_data!\n";
		for( std::size_t i = 0; i < n_data_fields(); ++i ){
			const std::string &s = data[i]->name;
			std::cerr << s << " ";
		}
		std::cerr << "\n";
		my_runtime_error( __FILE__, __LINE__,
		                  "Named data already in block_data" );
	}
	// Now you need to copy the data.
	data_field *cp = copy( &data_f );
	int index = data.size();
	data.push_back( cp );

	// Ignore some keys that are not unique for example:
	if( !is_legal_special_field( special_field ) ){
		return;
	}

	my_assert( __FILE__, __LINE__,
	           special_fields_by_index[special_field] == -1,
	           "Special field already set!" );

	special_fields_by_name[special_field]  = data_f.name;
	special_fields_by_index[special_field] = index;
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


void block_data::copy_meta( const block_data &o )
{
	tstep = o.tstep;
	N = o.N;
	N_types = o.N_types;
	atom_style = o.atom_style;
	dom = o.dom;
	top = o.top;

	for( int spec_field = block_data::ID;
	     spec_field < N_SPECIAL_FIELDS; ++spec_field ){
		const data_field *df = o.get_special_field( spec_field );
		if( !df ) continue;

		for( std::size_t i = 0; i < o.n_data_fields(); ++i ){
			const std::string name = o[i].name;
			if( name == df->name ){
				special_fields_by_name[spec_field]  = df->name;
				special_fields_by_index[spec_field] = i;
				break;
			}
		}
	}

}

void block_data::set_natoms( std::size_t new_size )
{
	N = new_size;
	for( data_field *d : data ){
		my_assert( __FILE__, __LINE__, d,
		           "Data field not initialised in attempted resize!" );
		d->resize(N);
	}
}


void block_data::set_special_field( const std::string &name, int field )
{
	int index = 0;
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

std::string block_data::get_special_field_name( int field ) const
{
	my_assert( __FILE__, __LINE__, is_legal_special_field( field ),
	           "Invalid field in get_special_field_name!" );

	return special_fields_by_name[field];
}




data_field *block_data::remove_field( const std::string &name,
                                      int &special_field )
{
	data_field *df = get_data_rw( name );
	if( !df ) return nullptr;

	// Delete the ptr from all vectors and maps.
	data.erase( std::find(data.begin(), data.end(), df) );
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
	}
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



// ******************   Non-member functions:    ************************
void swap( block_data &f, block_data &s )
{
	using std::swap;

	f.tstep = s.tstep;
	f.N = s.N;
	f.atom_style = s.atom_style;
	swap( f.dom, s.dom );
	swap( f.top, s.top );
	swap( f.data, s.data );
	swap( f.special_fields_by_name, s.special_fields_by_name );
	swap( f.special_fields_by_index, s.special_fields_by_index );
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


const data_field &block_data::operator[]( int i ) const
{
	return *data[i];
}

data_field &block_data::operator[]( int i )
{
	return *data[i];
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

} // namespace lammps_tools
