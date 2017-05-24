#include "dump_reader_lammps.hpp"

#include "dump_reader_lammps_plain.hpp"
#include "dump_reader_lammps_gzip.hpp"
#include "dump_reader_lammps_bin.hpp"


using namespace lammps_tools;
// using namespace readers;

namespace lammps_tools {

namespace readers {

bool is_int_data_field( const std::string &header )
{
	std::vector<std::string> ints = { "id", "mol", "type", "proc",
	                                  "ix", "iy", "iz" };
	return std::find( ints.begin(), ints.end(), header ) != ints.end();
}


void dump_reader_lammps::set_column_headers(
	const std::vector<std::string> &headers )
{
	column_headers = headers;
	column_header_types.resize( headers.size(), default_col_type );

	for( std::size_t i = 0; i < column_headers.size(); ++i ){
		// Assume that everything is a double, except for the cols
		// with is_int_data_field( col ) == true.
		if( is_int_data_field( column_headers[i] ) ){
			column_header_types[i] = data_field::INT;
		}else{
			column_header_types[i] = data_field::DOUBLE;
		}
	}
}

void dump_reader_lammps::set_column_header( std::size_t idx,
                                            const std::string &header )
{
	if( idx >= column_headers.size() ){
		column_headers.resize(idx+1);
		column_header_types.resize(idx+1);
	}
	my_assert( __FILE__, __LINE__, idx < column_headers.size(),
	           "Invalid index in set_column_header!" );
	column_headers[idx] = header;
	if( is_int_data_field( header ) ){
		column_header_types[idx] = data_field::INT;
	}
}

const std::vector<std::string> &dump_reader_lammps::get_column_headers() const
{
	return column_headers;
}


bool dump_reader_lammps::set_column_header_as_special( const std::string &header,
                                                       int special_field_type )
{
	my_assert( __FILE__, __LINE__,
	           is_legal_special_field( special_field_type ),
	           "Invalid special_field_type!" );

	if( std::find( column_headers.begin(), column_headers.end(),
	               header ) == column_headers.end() ){
		std::cerr << "Could not set " << header << " to type "
		          << special_field_type << ".\n";
		return false;
	}

	header_to_special_field[header] = special_field_type;
	if( is_special_field_int( special_field_type ) ){
		set_column_type( header, data_field::INT );
	}


	return true;
}


void dump_reader_lammps::add_custom_data_fields( std::vector<data_field*> &dfs,
                                                 block_data &b )
{
	for( data_field *df : dfs ){
		if( header_to_special_field.count( df->name ) ){
			int special_field_type =
				header_to_special_field[ df->name ];
			b.add_field( *df, special_field_type );
		}else{
			b.add_field( *df );
		}
		delete df;
	}
}



dump_reader_lammps *make_dump_reader_lammps( const std::string &fname,
                                             int fformat,
                                             std::vector<std::string> headers,
                                             int dump_style )
{
	dump_reader_lammps *reader = nullptr;
	if( fformat == FILE_FORMAT_PLAIN ){
		reader = new dump_reader_lammps_plain( fname, dump_style );
	}else if( fformat == FILE_FORMAT_BIN ){
		reader = new dump_reader_lammps_bin( fname, dump_style );
	}else if( fformat == FILE_FORMAT_GZIP ){
		reader = new dump_reader_lammps_gzip( fname, dump_style );
	}
	if( reader ) reader->set_column_headers( headers );
	return reader;
}

dump_reader_lammps *make_dump_reader_lammps( std::istream &input,
                                             std::vector<std::string> headers,
                                             int dump_style  )
{
	dump_reader_lammps *reader = nullptr;
	reader = new dump_reader_lammps_plain( input, dump_style );
	if( reader ) reader->set_column_headers( headers );
	return reader;

}

dump_reader_lammps *make_dump_reader_lammps( std::istream &input,
                                             int dump_style  )
{
	std::vector<std::string> empty;
	return make_dump_reader_lammps( input, empty, dump_style );
}

dump_reader_lammps *make_dump_reader_lammps( const std::string &fname, int fformat,
                                             int dump_style  )
{
	std::vector<std::string> empty;
	return make_dump_reader_lammps( fname, fformat, empty, dump_style );
}

void dump_reader_lammps::set_column_type( const std::string &header, int type )
{
	my_assert( __FILE__, __LINE__,
	           (type == data_field::INT || type == data_field::DOUBLE),
	           "Incorrect data type passed to set_column_type!");

	for( std::size_t i = 0; i < column_headers.size(); ++i ){
		if( column_headers[i] == header ){
			column_header_types[i] = type;
			return;
		}
	}

	std::stringstream ss;
	ss << "Header " << header << " not found!";
	my_warning( __FILE__, __LINE__, ss.str() );
}

int dump_reader_lammps::get_column_type( const std::string &header ) const
{
	for( std::size_t i = 0; i < column_headers.size(); ++i ){
		if( column_headers[i] == header ){
			return column_header_types[i];
		}
	}

	std::stringstream ss;
	ss << "Header " << header << " not found!";
	my_runtime_error( __FILE__, __LINE__, ss.str() );
	return -1;
}

void dump_reader_lammps::set_default_column_type( int type )
{
	if( type < data_field::DOUBLE || type > data_field::INT ){
		std::cerr << "Ignoring unrecognized type " << pretty_type(type)
		          << ", default remains "
		          << pretty_type(default_col_type) << ".\n";
		return;
	}
	std::cerr << "Putting default column type to "
	          << pretty_type(type) << ".\n";
	default_col_type = type;
}


} // namespace readers

} // namespace lammps_tools
