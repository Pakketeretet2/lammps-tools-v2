#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lt_dump_reader.h"
#include "../../../cpp_lib/enums.hpp"

PYBIND11_PLUGIN(dump_reader_) {
	pybind11::module m("dump_reader_", "Exposes dump_reader through pybind11");

	pybind11::class_<lt_dump_reader_handle>(m, "dump_reader_handle")
		.def(pybind11::init<>())
		.def_readonly("dr",&lt_dump_reader_handle::dr)
		.def_readonly("file_format",&lt_dump_reader_handle::fformat)
		.def_readonly("dump_format",&lt_dump_reader_handle::dformat);

	m.def("new",       &lt_new_dump_reader);
	m.def("new_local", &lt_new_dump_reader_local);

	m.def("delete", &lt_delete_dump_reader);
	m.def("status", &lt_dump_reader_status);
	m.def("get_next_block", &lt_get_next_block);
	m.def("number_of_blocks", &lt_number_of_blocks);
	m.def("set_column_header", &lt_set_col_header);
	m.def("set_special_column", &lt_set_column_header_as_special);

	using namespace lammps_tools;

	pybind11::enum_<LT_DUMP_READER_STATUS>(m, "DUMP_READER_STATUS")
		.value("IS_GOOD", LT_DUMP_READER_STATUS::IS_GOOD)
		.value("AT_EOF",  LT_DUMP_READER_STATUS::AT_EOF)
		.value("IS_BAD",  LT_DUMP_READER_STATUS::IS_BAD)
		.value("POINTER_NULL", LT_DUMP_READER_STATUS::POINTER_NULL);

	pybind11::enum_<LT_DUMP_READER_FILE_FORMATS>(m, "FILE_FORMATS")
		.value("PLAIN", FILE_FORMAT_PLAIN)
		.value("GZIP",  FILE_FORMAT_GZIP)
		.value("BIN",   FILE_FORMAT_BIN)
		.value("UNSET", FILE_FORMAT_UNSET);

	pybind11::enum_<LT_DUMP_READER_DUMP_FORMATS>(m, "DUMP_FORMATS")
		.value("LAMMPS", DUMP_FORMAT_LAMMPS)
		.value("HOOMD",  DUMP_FORMAT_HOOMD)
		.value("NAMD",   DUMP_FORMAT_NAMD)
		.value("UNSET",  DUMP_FORMAT_UNSET);





	return m.ptr();
}
