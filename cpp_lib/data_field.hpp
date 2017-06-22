#ifndef DATA_FIELD_HPP
#define DATA_FIELD_HPP

/**
   \file data_field.hpp

   Definition of data_field and associated functions.
*/

#include <string>
#include <typeinfo>
#include <vector>

#include "my_assert.hpp"
#include "util.hpp"

namespace lammps_tools {


/**
   \brief Contains any type of data that could be in a
   dump file together with some human readable identifier.
*/
struct data_field
{
	/// Defines possible data types contained in data_field.
	enum types { DOUBLE = 0, ///< data_field contains doubles
	             INT         ///< data_field contains ints
	};

	/// Constructor that only sets the name of the data_field
	explicit data_field( const std::string &n ) : name(n) {}
	/// Empty destructor
	virtual ~data_field() {}

	/// Returns the type of the underlying data. See \p types.
	virtual int type()     const = 0;
	/// Returns the size of the data_field
	virtual std::size_t size() const = 0;
	/// Resizes so that the data_field can accomodate N values.
	virtual void resize( std::size_t N ) = 0;

	/**
	   Friend swap function.

	   \param f The first block_data to swap with second.
	   \param s The second block_data to swap with first.
	*/
	friend void swap( data_field &f, data_field &s );

	std::string name; ///< The name of the header.
};


const char *pretty_type( int type );



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
	explicit data_field_der( const std::string &n )
		: data_field( n ), data()
	{}

	/**
	   Constructor that takes a name and a size.

	   Sets up name and size of data.
	*/
	data_field_der( const std::string &n, std::size_t size )
		: data_field( n ), data( size )
	{ }

	/**
	   Copy constructor from pointer.

	   Sets up name and complete data.
	*/
	explicit data_field_der( const data_field_der<T, TYPE> *other )
		: data_field( other->name ), data( other->size() )
	{
		construct_from_other( *other );
	}

	/**
	   Copy constructor from reference.

	   Sets up name and complete data.
	*/
	explicit data_field_der( const data_field_der<T, TYPE> &other )
		: data_field( other.name ), data( other.size() )
	{
		std::cerr << "Attempting to construct through ref.\n";
		construct_from_other( other );
	}

	/**
	   Explicitly set the data from given vector.
	*/
	explicit data_field_der( const std::string &n,
	                         const std::vector<T> &vec )
		: data_field(n), data(vec)
	{
		// std::cerr << "Called copy-from-vector constructor.\n";
	}

	/**
	   Move vector.
	*/
	explicit data_field_der( const std::string &n,
	                         std::vector<T> &&vec )
		: data_field(n), data(vec)
	{
		// std::cerr << "Called move-from-vector constructor.\n";
	}



	/**
	   \brief Helper to ease construction from other, either ref or pointer
	*/
	void construct_from_other( const data_field_der<T, TYPE> &other )
	{
		my_assert( __FILE__, __LINE__, other.type() == TYPE,
		           "Type mismatch in constructor!" );
		/*
		const std::vector<T> &d_o = other.get_data();
		for( std::size_t i = 0; i < size(); ++i ){
			data[i] = d_o[i];
		}
		*/
		data = other.get_data();
	}

	/**
	   Assingment operator
	*/
	data_field_der &operator=( data_field_der<T, TYPE> other )
	{
		swap( *this, other );
		return *this;
	}



	/**
	   Empty destructor.
	*/
	virtual ~data_field_der(){}

	/**
	   Returns type indicator.

	   \returns type of data.
	*/
	virtual int type() const { return TYPE; }

	/**
	   Returns the size of the data.

	  \returns size of data.
	*/
	virtual std::size_t size() const { return data.size(); }

	/**
	   Resizes the size of the data.

	   \param N New size data should have.
	*/
	virtual void resize( std::size_t N ){ data.resize( N ); }

	/**
	   Provides mutable access to underlying data vector.
	*/
	T& operator[]( std::size_t idx ) { return data[idx]; }

	/**
	   Provides copy from data vector.
	*/
	const T& operator[]( std::size_t idx ) const { return data[idx]; }


	/**
	   Provides immutable access to underlying data vector.
	*/
	virtual const std::vector<T> &get_data() const { return data; }

	/**
	   Provides read/write access to underlying data vector.
	*/
	virtual std::vector<T> &get_data_rw() { return data; }

	/**
	   Provides immutable access to underlying data array
	*/
	virtual const T* get_data_ptr() const { return data.data(); }


	/**
	   Actual implementation of swap for each type.
	*/
	// template<typename TT, int TTYPE>
	friend void swap( data_field_der &f, data_field_der &s )
	{
		using std::swap;

		swap( f.name, s.name );
		swap( f.data, s.data );
	}


private:
	std::vector<T> data; ///< The data vector.
};



/// \typedef Data field specialised for doubles
typedef data_field_der<double, data_field::DOUBLE> data_field_double;

/// \typedef Data field specialised for ints
typedef data_field_der<int, data_field::INT>    data_field_int;


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
			my_logic_error( __FILE__, __LINE__,
			                "Unkown data field type encountered!" );
			return nullptr;
		}
		case data_field::INT:{
			const data_field_int *dn =
				static_cast<const data_field_int*>( d );
			my_assert( __FILE__, __LINE__, dn != nullptr,
			           "Static cast to int failed!" );
			return new data_field_int( dn );
		}
		case data_field::DOUBLE:{
			const data_field_double *dn =
				static_cast<const data_field_double*>( d );
			my_assert( __FILE__, __LINE__, dn != nullptr,
			           "Static cast to double failed!" );
			return new data_field_double( dn );
		}
	}
}



template <typename T> inline
const std::vector<T> &data_as( const data_field *df );

template <typename T> inline
std::vector<T> &data_as_rw( data_field *df );




template <> inline
std::vector<double> &data_as_rw<double>( data_field *df )
{
	using dfd = data_field_double;
	int type = df->type();
	my_assert( __FILE__, __LINE__, type == data_field::DOUBLE,
	           "Type mismatch in cast in data_as<double>!" );

	dfd *df_c = static_cast<dfd*>( df );
	return df_c->get_data_rw();
}

template <> inline
std::vector<int> &data_as_rw<int>( data_field *df )
{
	using dfi = data_field_int;
	int type = df->type();
	my_assert( __FILE__, __LINE__, type == data_field::INT,
	           "Type mismatch in cast in data_as<int>!" );

	dfi *df_c = static_cast<dfi*>( df );
	return df_c->get_data_rw();
}


template <> inline
const std::vector<double> &data_as<double>( const data_field *df )
{
	using dfd = data_field_double;
	int type = df->type();
	my_assert( __FILE__, __LINE__, type == data_field::DOUBLE,
	           "Type mismatch in cast in data_as<double>!" );

	const dfd *df_c = static_cast<const dfd*>( df );
	return df_c->get_data();
}

template <> inline
const std::vector<int> &data_as<int>( const data_field *df )
{
	using dfi = data_field_int;
	int type = df->type();
	my_assert( __FILE__, __LINE__, type == data_field::INT,
	           "Type mismatch in cast in data_as<int>!" );

	const dfi *df_c = static_cast<const dfi*>( df );
	return df_c->get_data();
}


/**
   \brief Creates a permutation that sorts the given data.

   \param df    The data field to generate the permutation for.
   \param comp  Comparison functor/function pointer/std::function

   \returns a vector containing the sorting permutation.
*/
template<typename T, typename comparator> inline
std::vector<std::size_t> get_data_permutation( const data_field *df,
                                               comparator comp )
{
	const std::vector<T> &vec = data_as<T>( df );
	return util::sort_permutation( vec, comp );
}



/**
   \brief Sorts the given data field according to given permutation vector.

   \param df  The data field to sort
   \param p   The permutation to sort along.
*/
template<typename T> inline
void sort_data_with_permutation( data_field *df,
                                 const std::vector<std::size_t> &p )
{
	std::vector<T> &vec = data_as_rw<T>( df );
	util::apply_permutation( vec, p );
}



void swap( data_field_double &f, data_field_double &s );
void swap( data_field_int &f, data_field_int &s );

} // namespace lammps_tools


#endif // DATA_FIELD_HPP
