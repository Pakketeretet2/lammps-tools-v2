#ifndef BLOCK_DATA_H
#define BLOCK_DATA_H

#include "domain.hpp"
#include "data_field.hpp"
#include "topology.hpp"
#include "types.hpp"

/**
   A class that represents information about a single time step.
*/
class block_data
{
public:
	/// Possible atom styles
	enum atom_styles { ATOMIC,
	                   MOLECULAR };
	
	bigint tstep;    ///< The current time step
	bigint N;        ///< The number of atoms
	int atom_style;  ///< The current atom style
	domain dom;      ///< Domain information
	topology top;    ///< Topology information

	/**
	   Returns a pointer to an underlying data field.

	   \param name The name of the data field.

	   \returns A pointer to the found data field, or nullptr if the
	            data field is not available.
	*/
	data_field *get_field( const std::string &name );
	
	/**
	   Returns a pointer to an underlying data field.

	   \param name The name of the data field.

	   \returns A non-negative integer that indicates the type of the
	            named data field, or -1 if the field could not be found.
	*/
	int get_field_type( const std::string &name );

	/**
	   Adds a data field to the data fields.

	   \param data The data field to add.
	*/
	void add_field( const data_field &data );

	/**
	   Simple, empty constructor.
	*/
	block_data();

	/**
	   Copy constructor.
	*/
	block_data( const block_data &o );
	
	/**
	   Destructor. Makes sure contents of data are properly deleted.
	*/
	~block_data();

	/**
	   Assignment operator performs deep copy.
	   Deliberately void to prevent assignment chaining.

	   \param o The block_data whose contents to copy.
	*/
	void operator=( const block_data &o );

	/**
	   Swaps contents of this with o's.

	   \param o The block_data to swap with.
	*/
	void swap( block_data &o ) throw(); // For copy-and-swap idiom.

	/**
	   Copies only the meta-data from block_data o.

	   \param o the block_data to copy from.
	*/
	void copy_meta( const block_data &o );

	/**
	   Resizes all data_fields.

	   \param N New size to give to data fields.
	*/
	void resize( std::size_t N );

	/**
	   Initialises all data_fields.

	   \param N New size to give to data fields.
	*/
	void init( std::size_t N );

	/**
	   Returns the names of the data fields contained
	*/
	std::vector<std::string> get_data_names() const;

	/**
	   Returns the number of data fields contained.
	*/
	std::size_t get_data_size() const;

private:
	/// A vector containing pointers to all data fields.
	std::vector<data_field*> data;
};

/*
block_data block_data_from_data_file( const char *fname, int &status );


void print_block_data_lmp( const block_data &b, std::ostream &o );
void print_block_data_lmp( const block_data &b, const std::string &s );
void print_block_data_lmp( const block_data &b, const std::string &s,
                           std::ios_base::openmode mode );
*/




#endif // BLOCK_DATA_H
