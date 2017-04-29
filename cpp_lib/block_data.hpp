#ifndef BLOCK_DATA_HPP
#define BLOCK_DATA_HPP

#include "atom_type_info.hpp"
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
	int    N_types;  ///< The number of atom types
	int atom_style;  ///< The current atom style
	
	domain dom;      ///< Domain information
	topology top;    ///< Topology information
	

	/**
	   Simple, empty constructor.
	*/
	block_data();

	/**
	   Constructor that sets the expected number of atoms
	*/
	block_data( std::size_t n_atoms );

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

	   \param o The block_data whose contents to copy.

	   \returns Ref to this.
	*/
	block_data &operator=( block_data o );

	/**
	   Swaps contents of f and s.

	   \param f The first block_data to swap with second.
	   \param s The second block_data to swap with first.
	*/
	friend void swap( block_data &f, block_data &s );

	/**
	   Copies only the meta-data from block_data o.

	   \param o the block_data to copy from.
	*/
	void copy_meta( const block_data &o );


	/**
	   Returns a pointer to an underlying data field.

	   \warning It is _your_ responsibility to make sure
	            the data is actually present (by calling
	            get_data or checking for nullptr) before
	            using the data!

	   \param name The name of the data field.

	   \returns A ptr to the data field, or nullptr
                    if the field is not found.
	*/
	data_field *get_data_rw( const std::string &name );
	
	/**
	   Returns a pointer to an underlying data field.

	   \warning It is _your_ responsibility to make sure
	            the data is actually present (by calling
	            get_data or checking for nullptr) before
	            using the data!

	   \param name The name of the data field.

	   \returns A ptr to the data field, or nullptr
                    if the field is not found.
	*/
	const data_field *get_data( const std::string &name ) const;
	
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
	   Resizes all data_fields.

	   \param N New size to give to data fields.
	*/
	void set_natoms( std::size_t N );

	/**
	   Returns the names of the data fields contained
	*/
	std::vector<std::string> get_data_names() const;

	/**
	   Returns the number of data fields contained.
	*/
	std::size_t n_data_fields() const;

private:
	/// A vector containing pointers to all data fields.
	std::vector<data_field*> data;
};





#endif // BLOCK_DATA_HPP
