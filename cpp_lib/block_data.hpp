#ifndef BLOCK_DATA_HPP
#define BLOCK_DATA_HPP

/**
   \file block_data.hpp

   Definitions of block_data type and associated functions.
*/


#include "atom_type_info.hpp"
#include "domain.hpp"
#include "data_field.hpp"
#include "topology.hpp"
#include "types.hpp"

#include <map>
#include <string>




namespace lammps_tools {

/**
   A class that represents information about a single time step.
*/
class block_data
{
public:
	/// Possible atom styles
	enum atom_styles { ATOMIC = 0, ///< Atomic data only
	                   MOLECULAR   ///< Data on molecules too
	};

	/// enumerates "special" fields.
	enum special_fields { UNKNOWN = -1, ///< Use as error flag
	                      ID = 0,       ///< Atom ids
	                      TYPE,         ///< Atom types
	                      MOL,          ///< Molecule id
	                      X, 	    ///< x coordinate
	                      Y, 	    ///< y coordinate
	                      Z, 	    ///< z coordinate
	                      VX,	    ///< x velocity
	                      VY,	    ///< y velocity
	                      VZ,	    ///< z velocity
	                      IX,	    ///< x image flag
	                      IY,	    ///< y image flag
	                      IZ	    ///< z image flag
	};

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
	   Returns a pointer to an underlying data field.

	   \param name The name of the data field.

	   \returns A non-negative integer that indicates the type of the
	            named data field, or -1 if the field could not be found.
	*/
	int get_field_type( const std::string &name );

	/**
	   Adds a data field to the data fields.

	   \note The data provided is copied into the block_data.

	   \param data The data field to add.
	*/
	void add_field( const data_field &data, int special_field = UNKNOWN );

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
	   Returns the names of the data fields contained
	*/
	std::vector<std::string> get_data_names() const;

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


private:
	/// A vector containing pointers to all data fields.
	std::vector<data_field*> data;

	/// Contains a mapping of "special" fields to their respective names
	std::map<int, std::string> special_fields_by_name;

	/// Contains a mapping of "special" fields to their respective indices
	std::map<int, int> special_fields_by_index;
};



/**
   \brief Copies commonly used data fields from block_data.

   \param b        the block_data to extract the fields from.
   \param fields   in order the names of the id, type, x, y and z field
   \param id       will contain the id data vector
   \param type     will contain the type data vector
   \param x        will contain the x data vector
   \param y        will contain the y data_vector
   \param z        will contain the z data_vector

   If any of the fields cannot be found, all vectors shall be unchanged.

   \returns true on success (all fields found), false otherwise.
*/
bool grab_common_fields( const block_data &b,
                         const std::vector<std::string> &fields,
                         std::vector<int> &id, std::vector<int> &type,
                         std::vector<double> &x, std::vector<double> &y,
                         std::vector<double> &z );

/**
   \brief Grab data field by name and store it in vector of given type.

   \param b      block data whose field to grab
   \param field  data field by name
   \param vec    vector to store data in.

   \returns true on success, false on failure.
*/
template <typename T> inline
bool grab_field_as( const block_data &b, const std::string &field,
                    std::vector<T> &vec )
{
	my_logic_error( __FILE__, __LINE__,
	                "grab_field_as only supports int and double!" );
}

/// Specialisation of grab_field_as for double
template <> inline
bool  grab_field_as<double>( const block_data &b, const std::string &field,
                             std::vector<double> &vec )
{
	const data_field *df = b.get_data( field );
	if( !df ) return false;

	my_assert( __FILE__, __LINE__, df->type() == data_field::DOUBLE,
	           "Incorrect type for grab_field_as<double>!" );
	vec = data_as<double>( df );
	return true;
}

/// Specialisation of grab_field_as for double
template <> inline
bool  grab_field_as<int>( const block_data &b, const std::string &field,
                          std::vector<int> &vec )
{
	const data_field *df = b.get_data( field );
	if( !df ) return false;

	my_assert( __FILE__, __LINE__, df->type() == data_field::INT,
	           "Incorrect type for grab_field_as<int>!" );
	vec = data_as<int>( df );
	return true;
}


} // namespace lammps_tools

#endif // BLOCK_DATA_HPP
