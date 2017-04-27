#ifndef DATA_FIELD_HPP
#define DATA_FIELD_HPP

#include <string>
#include <vector>

#include "my_assert.hpp"

/**
   This is basically any type of data that could be in a
   dump file together with some human readable identifier.
*/
struct data_field
{
	enum types { DOUBLE = 0,
	             INT };
	
	data_field( const std::string &n ) : name(n) {}
	virtual ~data_field(){}
	std::string name;

	virtual int type()     const = 0;
	virtual std::size_t size() const = 0;
};



/**
   Derived classes actually provide acceptable data types.

   \param T     Type of the data contained.
   \param TYPE  An int that represents the type used (see data_field::types).

   \inherits data_field
*/
template <typename T, int TYPE>
struct data_field_der : public data_field
{
	/**
	   Constructor that takes a name.
	   
	   Sets up name only.
	*/
	data_field_der( const std::string &n ) : data_field( n ) {}

	/**
	   Constructor that takes a name and a size.
	   
	   Sets up name and size of data.
	*/
	data_field_der( const std::string &n, std::size_t size ) : data_field( n )
	{
		data.resize(size);
	}

	/**
	   Copy constructor from pointer.
	   
	   Sets up name and complete data.
	*/
	data_field_der( const data_field_der<T, TYPE> *other ) : data_field( other->name )
	{
		my_assert( other->type() == TYPE, "Type mismatch in constructor!",
		           __FILE__, __LINE__ );
		data.resize( other->size() );
		const std::vector<T> &d_o = other->get_data();
		for( std::size_t i = 0; i < size(); ++i ){
			data[i] = d_o[i];
		}
	}
	
	/**
	   Empty destructor.
	*/
	virtual ~data_field_der(){}

	/**
	   Returns type indicator.
	*/
	virtual int type() const { return TYPE; }

	/**
	   Returns the size of the data.
	*/
	virtual std::size_t size() const { return data.size(); }

	/**
	   Provides mutable access to underlying data vector.
	*/
	T& operator[]( std::size_t idx ) { return data[idx]; }

	/**
	   Provides immutable access to underlying data vector.
	*/
	virtual const std::vector<T> &get_data() const { return data; }
	
private:
	std::vector<T> data; ///< The data vector.
};


/// Data field containing doubles.
typedef data_field_der<double, data_field::DOUBLE> data_field_double;

/// Data field containing integers.
typedef data_field_der<int,    data_field::INT>    data_field_int;


/**
   Creates a copy of d, properly respecting its derived type.

   \param d The data to copy.
   
   \returns A pointer to a heap-allocated data_field.
*/
inline data_field *copy( const data_field *d )
{
	int type = d->type();
	switch( type ){
		default:{
			logic_error( "Unkown data field type encountered!",
			             __FILE__, __LINE__ );
			return nullptr;
		}
		case data_field::INT:{
			const data_field_int *dn =
				dynamic_cast<const data_field_int*>( d );
			my_assert( dn != nullptr, "Dynamic cast to int failed!",
			           __FILE__, __LINE__ );
			return new data_field_int( dn );
		}
		case data_field::DOUBLE:{
			const data_field_double *dn =
				dynamic_cast<const data_field_double*>( d );
			my_assert( dn != nullptr, "Dynamic cast to double failed!",
			           __FILE__, __LINE__ );
			return new data_field_double( dn );
		}
	}
}



#endif // DATA_FIELD_HPP
