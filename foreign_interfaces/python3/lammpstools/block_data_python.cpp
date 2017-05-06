#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "lammps_tools.h"

PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);


PYBIND11_PLUGIN(block_data) {
	pybind11::module m("block_data", "Exposes lt_block_data through pybind11");

	
	
	pybind11::class_<lt_block_data_handle>(m, "block_data_handle")
		.def(pybind11::init<>())
		.def("time_step"   , &lt_block_data_handle::time_step)
		.def("n_atoms"     , &lt_block_data_handle::n_atoms)
		.def("n_atom_types", &lt_block_data_handle::n_atom_types)
		.def("atom_style"  , &lt_block_data_handle::atom_style);

	pybind11::class_<lammps_tools::data_field_double>(m, "data_field_double")
		.def(pybind11::init<const std::string &>());

	m.def("get_data", &lt_get_data,
	      "Gets named data from block_data_handle");
	m.def("data_as_double", &lt_data_as_double,
	      "Gets named data from block_data_handle interpreted as doubles");

	
	return m.ptr();
}

