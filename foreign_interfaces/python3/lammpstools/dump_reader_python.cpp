#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lt_dump_reader.h"

PYBIND11_PLUGIN(dump_reader_) {
	pybind11::module m("dump_reader_", "Exposes dump_reader through pybind11");

	pybind11::class_<lt_dump_reader_handle>(m, "dump_reader_handle")
		.def(pybind11::init<>())
		.def_readonly("dr",&lt_dump_reader_handle::dr)
		.def_readonly("file_format",&lt_dump_reader_handle::fformat)
		.def_readonly("dump_format",&lt_dump_reader_handle::dformat);

	m.def("new_dump_reader", &lt_new_dump_reader);
	m.def("delete_dump_reader", &lt_delete_dump_reader);
	m.def("dump_reader_status", &lt_dump_reader_status);
	m.def("get_next_block", &lt_get_next_block);
	m.def("number_of_blocks", &lt_number_of_blocks);
	m.def("set_column_header", &lt_set_col_header);

	pybind11::enum_<LT_DUMP_READER_STATUS>(m, "DUMP_READER_STATUS")
		.value("IS_GOOD", LT_DUMP_READER_STATUS::IS_GOOD)
		.value("AT_EOF",  LT_DUMP_READER_STATUS::AT_EOF)
		.value("IS_BAD",  LT_DUMP_READER_STATUS::IS_BAD)
		.value("POINTER_NULL", LT_DUMP_READER_STATUS::POINTER_NULL);




	return m.ptr();
}
