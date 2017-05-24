#ifndef BLOCK_DATA_HPP
#define BLOCK_DATA_HPP

/**
   \file block_data.hpp

   Definitions of block_data type and associated functions.
*/


#include "atom_type_info.hpp"
#include "domain.hpp"
#include "data_field.hpp"
#include "enums.hpp"
#include "topology.hpp"
#include "types.hpp"

#include <map>
#include <string>




namespace lammps_tools {


/**
   Returns true if the special_field should be treated as int instead of double
*/
bool is_special_field_int( int special_field );


/**
   A class that represents information about a single time step.

   \warning This class has a custom copy operator, so whenever
            new members are added, make sure they are properly
            included in the copy operation as well.
*/
class block_data
{
public:
	/// enumerates "special" fields.
	enum special_fields { UNKNOWN = -1,     ///< Use as error flag
	                      ID = 0,           ///< Atom ids
	                      TYPE,             ///< Atom types
	                      MOL,              ///< Molecule id
	                      X, 	        ///< x coordinate
	                      Y, 	        ///< y coordinate
	                      Z, 	        ///< z coordinate
	                      VX,	        ///< x velocity
	                      VY,	        ///< y velocity
	                      VZ,	        ///< z velocity
	                      IX,	        ///< x image flag
	                      IY,	        ///< y image flag
	                      IZ,	        ///< z image flag

	                      /// Fake entry, counts number of special fields.
	                      N_SPECIAL_FIELDS
	};

	// Public members:
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
	explicit block_data( std::size_t n_atoms );

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

	   \warning It is _your_ responsibility to make sure the data
	            is actually present (by calling get_field_type or
		    checking for nullptr) before using the data!

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
	   Adds a data field to the data fields.

	   \note The data provided is copied into the block_data.

	   \param data The data field to add.
	   \parma special_field_type the type of specialty of data.
	*/
	void add_field( const data_field &data,
	                int special_field_type = block_data::UNKNOWN );

	/**
	   \brief Removes the named data field from the block_data.

	   This returns a pointer to the field so that it can be cleaned up.

	   \param[in]  name           the name of the data field to remove
	   \param[out] special_field  the special_field type of data

	   \returns a shared pointer to the data field, or nullptr if not found
	*/
	data_field *remove_field( const std::string &name, int &special_field );

	/**
	   Resizes all data_fields.

	   \param N New size to give to data fields.
	*/
	void set_natoms( std::size_t N );

	/**
	   Returns the number of data fields contained.
	*/
	std::size_t n_data_fields() const;

	/**
	   Sets the given name as given special field.
	*/
	void set_special_field( const std::string &name, int field );

	/**
	   Get name of given special field, or empty string if it is not set.

	   \param field  Identifier of the special field (see special_fields)

	   \return The name of the data_field that is special, or
                   empty string if the special field is not set.
	*/
	std::string get_special_field_name( int field ) const;

	/// Get read/write pointer to special data field of given kind.
	data_field *get_special_field_rw( int field );

	/// Get read-only pointer to special data field of given kind.
	const data_field *get_special_field( int field ) const;

	/// Sorts the data along given header with given comparator.
	template <typename comparator>
	void sort_along( const std::string &header, const comparator &comp );

	/// Sorts the data along given header with standard '<' comparator.
	void sort_along( const std::string &header );

	/// Grab fields by index:
	const data_field &operator[]( int i ) const;

	/// Grab fields by index:
	data_field &operator[]( int i );

private:
	/// A vector containing pointers to all data fields.
	std::vector<data_field*> data;

	/// Contains a mapping of "special" fields to their respective names
	std::vector<std::string> special_fields_by_name;

	/// Contains a mapping of "special" fields to their respective indices
	std::vector<int> special_fields_by_index;
};




inline
void block_data::sort_along( const std::string &header )
{
	auto comp = []( double x, double y ) { return x < y; };
	sort_along( header, comp );
}

template <typename comparator> inline
void block_data::sort_along( const std::string &header, const comparator &comp )
{
	// First find the sorting permutation, then apply it to all fields.
	const data_field *df = get_data( header );
	my_assert( __FILE__, __LINE__, df,
	           "Data field for sort does not exist!" );

	std::vector<std::size_t> p;
	if( df->type() == data_field::DOUBLE ){
		p = get_data_permutation<double>( df, comp );
	}else if( df->type() == data_field::INT ){
		p = get_data_permutation<int>( df, comp );
	}else{
		my_logic_error( __FILE__, __LINE__,
		                "Unkown data type in block_data::sort_along!" );
	}

	for( data_field *df : data ){
		if( df->type() == data_field::DOUBLE ){
			sort_data_with_permutation<double>( df, p );
		}else if( df->type() == data_field::INT ){
			sort_data_with_permutation<int>( df, p );
		}
	}
}


/**
   \brief Checks whether special_field is _legal_, meaning it is an
   understood special_field or unknown

   \param special_field the special_field to test.

   \returns true if special_field is legal, false otherwise.
*/
bool is_legal_special_field( int special_field );




} // namespace lammps_tools


#endif // BLOCK_DATA_HPP
