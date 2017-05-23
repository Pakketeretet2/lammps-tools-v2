#ifndef BLOCK_DATA_ACCESS_HPP
#define BLOCK_DATA_ACCESS_HPP

/**
   \file block_data_access.hpp

   Provides utility accessors for grabbing data from a block_data object.
*/


#include "block_data.hpp"

namespace lammps_tools {

// ******************* non-member functions:  *******************************

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
	return false;
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


// Convenience functions:
inline const std::vector<double> &get_x( const block_data &b )
{
	return data_as<double>( b.get_special_field( block_data::X ) );
}

inline const std::vector<double> &get_y( const block_data &b )
{
	return data_as<double>( b.get_special_field( block_data::Y ) );
}

inline const std::vector<double> &get_z( const block_data &b )
{
	return data_as<double>( b.get_special_field( block_data::Z ) );
}

inline const std::vector<int> &get_id( const block_data &b )
{
	return data_as<int>( b.get_special_field( block_data::ID ) );
}

inline const std::vector<int> &get_type( const block_data &b )
{
	return data_as<int>( b.get_special_field( block_data::TYPE ) );
}

inline const std::vector<int> &get_mol( const block_data &b )
{
	return data_as<int>( b.get_special_field( block_data::MOL ) );
}




inline std::vector<double> &get_x_rw(  block_data &b )
{
	return data_as_rw<double>( b.get_special_field_rw( block_data::X ) );
}

inline std::vector<double> &get_y_rw(  block_data &b )
{
	return data_as_rw<double>( b.get_special_field_rw( block_data::Y ) );
}

inline std::vector<double> &get_z_rw(  block_data &b )
{
	return data_as_rw<double>( b.get_special_field_rw( block_data::Z ) );
}

inline std::vector<int> &get_id_rw(  block_data &b )
{
	return data_as_rw<int>( b.get_special_field_rw( block_data::ID ) );
}

inline std::vector<int> &get_type_rw(  block_data &b )
{
	return data_as_rw<int>( b.get_special_field_rw( block_data::TYPE ) );
}

inline std::vector<int> &get_mol_rw(  block_data &b )
{
	return data_as_rw<int>( b.get_special_field_rw( block_data::MOL ) );
}


} // namespace lammps_tools

#endif // BLOCK_DATA_ACCESS_HPP
