#include "dump_reader_lammps.hpp"

namespace dump_readers {

bool is_int_data_field( const std::string header )
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

const std::vector<std::string> &dump_reader_lammps::get_column_headers() const
{
	return column_headers;
}




} // namespace dump_readers
