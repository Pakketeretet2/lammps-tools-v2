#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lt_block_data.h"
#include "lt_data_field.h"


PYBIND11_PLUGIN(data_field_) {

	pybind11::module m("data_field_", "Exposes data_field through pybind11");

	pybind11::class_<lt_data_field_handle>(m, "data_field_handle")
		.def(pybind11::init<>());

	m.def("get_size", &lt_data_field_size, "Returns the size of data field.");
	m.def("get_type", &lt_data_type, "Returns the type of data field.");
	m.def("get_name", &lt_data_name, "Returns the name of data field.");
	m.def("as_float", &lt_data_as_double_vec, "Returns the data as floats.");
	m.def("as_int",   &lt_data_as_int_vec,    "Returns the data as ints.");

	return m.ptr();
}
