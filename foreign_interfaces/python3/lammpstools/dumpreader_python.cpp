#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lammps_tools.h"

#include <cstring>

PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);

void lt_set_col_headers_helper( lt_dump_reader_handle drh, int n,
                                const std::vector<std::string> &headers )
{
	std::size_t idx = 0;
	for( const std::string &s : headers ){
		lt_set_col_header( drh, idx, s.c_str() );
		++idx;
	}

}

void list_of_strings( const std::vector<std::string> &list )
{
	for( const std::string &s : list ){
		std::cerr << s << "\n";
	}
}


PYBIND11_PLUGIN(dumpreader) {
	pybind11::module m("dump_reader", "Exposes lt_dump_reader through pybind11");
	
	pybind11::class_<lt_dump_reader_handle>(m, "dump_reader_handle");

	// Dump reader stuff:
	m.def("new_dump_reader", &lt_new_dump_reader,
	      pybind11::return_value_policy::copy,
	      "Initialises a new dump_reader and returns a handle to it.");
	m.def("delete_dump_reader", &lt_delete_dump_reader,
	      "Deletes given dump_reader.");
	m.def("get_next_block", &lt_get_next_block,
	      "Reads in a new block to given block_data_handle.");
	m.def("set_column_headers", &lt_set_col_headers_helper,
	      "Sets the column headers to given values.");
	m.def("print_list_of_strings", &list_of_strings,
	      "Prints a given list of strings.");

	return m.ptr();
}

