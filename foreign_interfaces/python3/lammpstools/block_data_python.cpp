#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lt_block_data.h"


PYBIND11_PLUGIN(block_data_) {
	pybind11::module m("block_data_", "Exposes block_data through pybind11");

	pybind11::class_<lt_block_data_handle>(m, "block_data_handle")
		.def(pybind11::init<>())
		.def(pybind11::init<lt_block_data_handle &>())
		.def("time_step", &lt_block_data_handle::time_step)
		.def("n_atoms", &lt_block_data_handle::n_atoms);

	// m.def("get_data", &lt_get_data, "Get block data from dump reader handle" );
	m.def("special_field_double", &lt_special_field_double,
	      "Get special field interpreted as double" );
	m.def("special_field_int",    &lt_special_field_int,
	      "Get special field interpreted as int" );


	return m.ptr();
}
