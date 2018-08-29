#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "lt_block_data.h"
#include "lt_data_field.h"

#include "make_vectors_opaque.hpp"

PYBIND11_PLUGIN(data_field_) {
	pybind11::module m("data_field_", "Exposes data_field through pybind11");


	pybind11::class_<lt_data_field_handle>(m, "data_field_handle")
		.def(pybind11::init<>());


	m.def("get_size", &lt_data_field_size, "Returns the size of data field.");
	m.def("get_type", &lt_data_field_type, "Returns the type of data field.");
	m.def("get_name", &lt_data_field_name, "Returns the name of data field.");
	m.def("as_float", &lt_data_as_double_vec, "Returns the data as floats.");
	m.def("as_int",   &lt_data_as_int_vec,    "Returns the data as ints.");

	// Expose some ways to modify the data from Python:
	m.def("set_size", &lt_data_field_set_size, "Sets the size of data field.");
	m.def("set_name", &lt_data_field_set_name, "Sets the name of data field.");
	m.def("new_data_field", &lt_new_data_field, "Creates a new data field.");
	m.def("delete_data_field", &lt_delete_data_field, "Deletes data field.");

	m.def("get_indexed_int_data", &lt_data_field_get_indexed_int_data,
	      "Grab single value as int by index." );
	m.def("get_indexed_double_data", &lt_data_field_get_indexed_double_data,
	      "Grab single value as double by index." );

	m.def("set_indexed_int_data", &lt_data_field_set_indexed_int_data,
	      "Set single value as int by index." );
	m.def("set_indexed_double_data", &lt_data_field_set_indexed_double_data,
	      "Set single value as double by index." );



	pybind11::enum_<LT_DATA_FIELD_TYPES>(m, "TYPES")
		.value("DOUBLE", LT_DATA_FIELD_TYPES::DATA_FIELD_DOUBLE)
		.value("INT",    LT_DATA_FIELD_TYPES::DATA_FIELD_INT);


	return m.ptr();
}
