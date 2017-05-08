#include "dump_reader_lammps.hpp"

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
}

void dump_reader_lammps::set_column_header( std::size_t idx,
                                            const std::string &header )
{
	if( idx >= column_headers.size() ){
		column_headers.resize(idx+1);
	}
	my_assert( __FILE__, __LINE__, idx < column_headers.size(),
	           "Invalid index in set_column_header!" );
	column_headers[idx] = header;
}

const std::vector<std::string> &dump_reader_lammps::get_column_headers() const
{
	return column_headers;
}


bool dump_reader_lammps::set_column_header_as_special( const std::string &header,
                                                       int special_field_type )
{
	my_assert( __FILE__, __LINE__,
	           special_field_type >= block_data::ID &&
	           special_field_type <= block_data::IZ,
	           "Invalid special_field_type!" );

	if( std::find( column_headers.begin(), column_headers.end(),
	               header ) == column_headers.end() ){
		std::cerr << "Could not set " << header << " to type "
		          << special_field_type << ".\n";
		return false;
	}

	header_to_special_field[header] = special_field_type;
	return true;
}


void dump_reader_lammps::add_custom_data_fields( std::vector<data_field*> &dfs,
                                                 block_data &b )
{
	for( data_field *df : dfs ){
		int special_field_type;
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

} // namespace readers

} // namespace lammps_tools
