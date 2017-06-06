#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "lt_block_data.h"

PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);

PYBIND11_PLUGIN(block_data_) {
	pybind11::module m("block_data_", "Exposes block_data through pybind11");

	pybind11::bind_vector<std::vector<int>>(m, "VectorInt");
	pybind11::bind_vector<std::vector<double>>(m, "VectorDouble");

	pybind11::class_<lammps_tools::block_data>(m, "block_data")
		.def(pybind11::init<>());

	pybind11::class_<lt_block_data_handle>(m, "block_data_handle")
		.def(pybind11::init<>())
		.def(pybind11::init<lt_block_data_handle &>())
		.def("time_step", &lt_block_data_handle::time_step)
		.def("n_atoms", &lt_block_data_handle::n_atoms)
		.def("get_const_ref", &lt_block_data_handle::get_const_ref);

	// m.def("get_data", &lt_get_data, "Get block data from dump reader handle" );
	m.def("has_special_field", &lt_has_special_field,
	      "Checks if block has special field" );
	m.def("special_field_double", &lt_special_field_double,
	      "Get special field interpreted as double" );
	m.def("special_field_int",    &lt_special_field_int,
	      "Get special field interpreted as int" );

	m.def("n_data_fields", &lt_n_data_fields,
	      "Get number of data fields.");
	m.def("data_by_index", &lt_data_by_index,
	      "Get nth data field.");
	m.def("data_by_name", &lt_data_by_name,
	      "Get data field by name");

	pybind11::enum_<lammps_tools::block_data::special_fields>(m, "SPECIAL_COLS")
		.value("ID", lammps_tools::block_data::ID)
		.value("TYPE", lammps_tools::block_data::TYPE)
		.value("MOL", lammps_tools::block_data::MOL)
		.value("X", lammps_tools::block_data::X)
		.value("Y", lammps_tools::block_data::Y)
		.value("Z", lammps_tools::block_data::Z)
		.value("VX", lammps_tools::block_data::VX)
		.value("VY", lammps_tools::block_data::VY)
		.value("VZ", lammps_tools::block_data::VZ)
		.value("IX", lammps_tools::block_data::IX)
		.value("IY", lammps_tools::block_data::IY)
		.value("IZ", lammps_tools::block_data::IZ);


	return m.ptr();
}
