#ifndef DUMP_READER_LAMMPS_HPP
#define DUMP_READER_LAMMPS_HPP

// LAMMPS dump readers come in three flavours: Plain, gzip and bin.
// They all share this general interface.

#include <iosfwd>

#include "dump_reader.hpp"

namespace lammps_tools {

class block_data;

namespace readers {

class dump_reader_lammps : public dump_reader
{
public:
	enum dump_styles { ATOMIC,
	                   CUSTOM };

	/// Empty constructor
	dump_reader_lammps(){}


	/// Empty destructor
	virtual ~dump_reader_lammps() {}

	/// Sets up a vector containing the expected columns.
	void set_column_headers( const std::vector<std::string> &headers );
	void set_column_header( std::size_t idx, const std::string &header );

	/// Returns a vector containing the expected column headers.
	const std::vector<std::string> &get_column_headers() const;

	/**
	   \brief Marks a given column header as a special field.
	   
	   \param header
	   \param special_field_type

	   If the header could not be found, block_data shall be unchanged.

	   \returns true if the header could be found, false otherwise.
	*/
	bool set_column_header_as_special( const std::string &header,
	                                   int special_field_type );


	/**
	   \brief Adds given vector of data fields to given block.
	   
	   \note This function deletes the data fields!

	   \param[in/out] dfs Data fields to add. They are deleted afterwards.
	   \param[out]    b   The block data to which the fields are added.
	*/
	void add_custom_data_fields( std::vector<data_field*> &dfs,
	                             block_data &b );
	
protected:
	std::map<std::string, int> header_to_special_field;
private:
	virtual int  get_next_block( lammps_tools::block_data &block ) = 0;
	virtual bool check_eof()  const = 0;
	virtual bool check_good() const = 0;

	std::vector<std::string> column_headers;
};


/// Checks if given column header corresponds to an integer quantity.
bool is_int_data_field( const std::string &header );

} // namespace readers

} // namespace lammps_tools

#endif // DUMP_READER_LAMMPS_HPP
