#include "block_data.hpp"
#include "my_assert.hpp"

block_data::block_data() : tstep( 0 ),
                           N( 0 ),
                           N_types( 1 ),
                           atom_style( ATOMIC )
{ }

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

void block_data::add_field( const data_field &data_f )
{
	my_assert( __FILE__, __LINE__, data_f.size() == N,
	           "Atom number mismatch on add_field! Call set_natoms first!");
	// Check if this name is already in block or not.
	my_assert( __FILE__, __LINE__, get_data( data_f.name ) == nullptr,
	           "Named data already in block_data" );

	
	// Now you need to copy the data.
	data_field *cp = copy( &data_f );
	data.push_back( cp );
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
