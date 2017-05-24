#ifndef DUMP_READER_LAMMPS_HPP
#define DUMP_READER_LAMMPS_HPP

// LAMMPS dump readers come in three flavours: Plain, gzip and bin.
// They all share this general interface.

#include <iosfwd>

#include "dump_reader.hpp"

namespace lammps_tools {

class block_data;

namespace readers {


// forward decl of implemented dump_reader_lammpses.
class dump_reader_lammps_plain;
class dump_reader_lammps_gzip;
class dump_reader_lammps_bin;


class dump_reader_lammps : public dump_reader
{
public:
	enum dump_styles { ATOMIC,
	                   CUSTOM,
	                   LOCAL };

	/// Empty constructor
	dump_reader_lammps( int dump_style )
		: dump_style(dump_style), default_col_type(data_field::DOUBLE),
		  header_to_special_field(), column_headers(),
		  column_header_types()
	{
		std::cerr << "Initiated dump_reader_lammps with "
		          << "dump_style " << dump_style << ".\n";
	}

	const int dump_style;

	/// Empty destructor
	virtual ~dump_reader_lammps() {}

	/// Sets up a vector containing the expected columns.
	void set_column_headers( const std::vector<std::string> &headers );

	void set_column_header( std::size_t idx, const std::string &header,
	                        int special_field_type = block_data::UNKNOWN );

	void set_column_type( const std::string &header, int type );
	void set_default_column_type( int type );

	/// Returns a vector containing the expected column headers.
	const std::vector<std::string> &get_column_headers() const;
	int get_column_type( const std::string &header ) const;

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


	/**
	   \brief Adds given vector of data fields to given block.

	   \note This function deletes the data fields!

	   \param[in/out] dfs Data fields to add. They are deleted afterwards.
	   \param[out]    b   The block data to which the fields are added.
	*/
	void add_local_data_fields( std::vector<data_field*> &dfs,
	                            block_data &b );

	int default_col_type; ///< Stores the default value assumed for columns.

protected:
	std::map<std::string, int> header_to_special_field;
private:
	virtual int  get_next_block( lammps_tools::block_data &block ) = 0;
	virtual bool check_eof()  const = 0;
	virtual bool check_good() const = 0;

	std::vector<std::string> column_headers;
	std::vector<int> column_header_types;

};


/// Checks if given column header corresponds to an integer quantity.
bool is_int_data_field( const std::string &header );





/**
   \brief Constructs a LAMMPS dump reader from file name.

   \param fname    Name of the dump file.
   \param fformat  File format (see \p FILE_FORMATS).
   \param header   Optionally you can pass the expected headers.
                   If the file is binary these are assumed to be
                   the column headers. For other file formats these
                   _have_ to match the column headers in the file
                   and are used for verification.

  \returns a pointer to the open dump_reader object, or nullptr on failure.

  \warning The dump_reader object is heap-allocated. Be sure to delete it
           or wrap it in a smart pointer.
*/
dump_reader_lammps *make_dump_reader_lammps( const std::string &fname, int fformat,
                                             std::vector<std::string> headers,
                                             int dump_style = dump_reader_lammps::CUSTOM );

/**
   \brief Constructs a LAMMPS dump reader from input file stream.

   \param fname    Name of the dump file.
   \param header   Optionally you can pass the expected headers.
                   If the file is binary these are assumed to be
                   the column headers. For other file formats these
                   _have_ to match the column headers in the file
                   and are used for verification.

  \returns a pointer to the open dump_reader object, or nullptr on failure.

  \warning The dump_reader object is heap-allocated. Be sure to delete it
           or wrap it in a smart pointer.
*/
dump_reader_lammps *make_dump_reader_lammps( std::istream &input,
                                             std::vector<std::string> headers,
                                             int dump_style = dump_reader_lammps::CUSTOM );


/**
   \brief Constructs a LAMMPS dump reader from file name.
   \overloads make_dump_reader_lammps
*/
dump_reader_lammps *make_dump_reader_lammps( const std::string &fname, int fformat,
                                             int dump_style = dump_reader_lammps::CUSTOM );


/**
   \brief Constructs a LAMMPS dump reader from input file stream.
   \overloads make_dump_reader_lammps
*/
dump_reader_lammps *make_dump_reader_lammps( std::istream &input,
                                             int dump_style = dump_reader_lammps::CUSTOM  );





} // namespace readers

} // namespace lammps_tools

#endif // DUMP_READER_LAMMPS_HPP
