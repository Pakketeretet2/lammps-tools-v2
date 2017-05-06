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

} // namespace readers

} // namespace lammps_tools
