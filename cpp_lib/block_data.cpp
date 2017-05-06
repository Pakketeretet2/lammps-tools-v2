#include "block_data.hpp"
#include "my_assert.hpp"

using namespace lammps_tools;

namespace lammps_tools {

block_data::block_data() : tstep( 0 ),
                           N( 0 ),
                           N_types( 1 ),
                           atom_style( ATOMIC )
{}

block_data::block_data( std::size_t n_atoms ) : tstep( 0 ),
                                                N( n_atoms ),
                                                N_types( 1 ),
                                                atom_style( ATOMIC )
{ }

block_data::block_data( const block_data &o )
	: tstep( o.tstep ),
	  N( o.N ),
	  N_types( o.N_types ),
	  atom_style( o.atom_style ),
	  data( o.n_data_fields() ),
	  dom( o.dom ),
	  top( o.top )

{
	my_assert( __FILE__, __LINE__, data.size() == o.n_data_fields(),
	           "Data size mismatch after copy!" );
	const std::vector<std::string> &names = o.get_data_names();
	for( std::size_t i = 0; i < o.n_data_fields(); ++i ){
		const std::string &name = names[i];
		data[i] = copy( o.get_data(name) );
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


int  block_data::get_field_type( const std::string &name )
{
	if( const data_field *df = get_data(name) ){
		return df->type();
	}
	return -1;
}

void block_data::add_field( const data_field &data_f, int special_field )
{
	my_assert( __FILE__, __LINE__, data_f.size() == N,
	           "Atom number mismatch on add_field! Call set_natoms first!");
	// Check if this name is already in block or not.
	my_assert( __FILE__, __LINE__, get_data( data_f.name ) == nullptr,
	           "Named data already in block_data" );

	// Now you need to copy the data.
	data_field *cp = copy( &data_f );
	int index = data.size();
	data.push_back( cp );

	// Ignore some keys that are not unique for example:
	if( special_field < ID || special_field > IZ ){
		return;
	}

	special_fields_by_name.insert(
		std::make_pair( special_field, data_f.name ) );
	special_fields_by_index.insert(
		std::make_pair( special_field, index ) );
}


block_data &block_data::operator=( block_data o )
{
	swap( *this, o );
	return *this;
}


std::vector<std::string> block_data::get_data_names() const
{
	std::vector<std::string> names( data.size() );
	std::size_t i = 0;
	for( const data_field *d : data ){
		names[i] = d->name;
		++i;
	}
	return names;
}

std::size_t block_data::n_data_fields() const
{
	return data.size();
}


void block_data::copy_meta( const block_data &o )
{
	tstep = o.tstep;
	N = o.N;
	atom_style = o.atom_style;
	dom = o.dom;
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
	my_assert( __FILE__, __LINE__, field >= ID && field <= IZ,
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
	my_assert( __FILE__, __LINE__, field >= ID && field <= IZ,
	           "Invalid field in get_special_field_name!" );
	auto name_it = special_fields_by_name.find( field );
	if( name_it != special_fields_by_name.end() ){
		return name_it->second;
	}else{
		return "";
	}
}




data_field *block_data::remove_field( const std::string &name,
                                      int &special_field )
{
	data_field *df = get_data_rw( name );
	if( !df ) return nullptr;

	// Delete the ptr from all vectors and maps.
	data.erase( std::find(data.begin(), data.end(), df) );
	bool found = false;
	for( auto const &ent : special_fields_by_name ){
		if( ent.second == name ){
			special_field = ent.first;
			found = true;
			break;
		}
	}
	if( found ){
		special_fields_by_name.erase( special_field );
		special_fields_by_index.erase( special_field );
		return df;
	}else{
		return nullptr;
	}
}

void swap( block_data &f, block_data &s )
{
	using std::swap;

	f.tstep = s.tstep;
	f.N = s.N;
	f.atom_style = s.atom_style;
	swap(f.dom, s.dom);
	swap(f.top, s.top);
	swap( f.data, s.data );

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


} // namespace lammps_tools
