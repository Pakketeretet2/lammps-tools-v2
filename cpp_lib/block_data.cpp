#include "block_data.hpp"
#include "my_assert.hpp"

block_data::block_data() : tstep( 0 ), N( 0 ), atom_style( ATOMIC ) {}
block_data::block_data( const block_data &o )
	: tstep( o.tstep ), N( o.N ), atom_style( o.atom_style )
{
	my_assert( o.N == data.size(), "Data size mismatch after copy!", __FILE__, __LINE__ );
	
}

block_data::~block_data()
{
	for( data_field *df : data ){
		delete df;
	}
}

data_field *block_data::get_field( const std::string &name )
{
	for( data_field *df : data ){
		if( df->name == name ){
			return df;
		}
	}
	return nullptr;
}

int  block_data::get_field_type( const std::string &name )
{
	if( data_field *df = get_field(name) ){
		return df->type();
	}
	return -1;
}

void block_data::add_field( const data_field &data_f )
{

	// Now you need to copy the data.
	
	
	//data.push_back( data_f );
}


void block_data::operator=( const block_data &o )
{
	tstep = o.tstep;
	N = o.N;
	atom_style = o.atom_style;
	dom = o.dom;
	top = o.top;

	// You need to copy the data.
	
}


std::vector<std::string> block_data::get_data_names() const
{
	std::vector<std::string> names( data.size() );
	std::size_t i = 0;
	for( const data_field *d : data ){
		names[i] = d->name;
	}
	return names;
}

std::size_t block_data::get_data_size() const
{
	return data.size();
}


/*
void block_data::swap( block_data &o ) throw();
void block_data::copy_meta( const block_data &o );
void block_data::resize( std::size_t N );
void block_data::init( std::size_t N );
*/
