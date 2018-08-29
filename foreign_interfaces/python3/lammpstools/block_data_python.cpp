#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "../../../c_interface/lt_block_data.h"

#include "make_vectors_opaque.hpp"


const void *vector_int_to_cptr( const std::vector<int> &v )
{
	return v.data();
}

const void *vector_double_to_cptr( const std::vector<double> &v )
{
	return v.data();
}



PYBIND11_PLUGIN(block_data_) {
	pybind11::module m("block_data_", "Exposes block_data through pybind11");

	pybind11::bind_vector<std::vector<int>>(m, "VectorInt");
	pybind11::bind_vector<std::vector<double>>(m, "VectorDouble");
	pybind11::bind_vector<std::vector<std::complex<double>>>(m, "VectorComplexDouble");
	pybind11::bind_vector<std::vector<std::complex<int>>>(m, "VectorComplexInt");

	pybind11::class_<lammps_tools::block_data>(m, "block_data")
		.def(pybind11::init<>());

	pybind11::class_<lt_block_data_handle>(m, "block_data_handle")
		.def(pybind11::init<>())
		.def(pybind11::init<lt_block_data_handle &>())
		.def("time_step", &lt_block_data_handle::time_step)
		.def("n_atoms", &lt_block_data_handle::n_atoms)
		.def("n_types", &lt_block_data_handle::n_types)
		.def("get_const_ref", &lt_block_data_handle::get_const_ref)
		.def("get_ptr", &lt_block_data_handle::get_ptr)
		.def("set_n_atoms", &lt_block_data_handle::set_n_atoms)
		.def("set_n_types", &lt_block_data_handle::set_n_types)
		.def("set_atom_style", &lt_block_data_handle::set_atom_style);

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

	// Allow grabbing the raw pointer of VectorInt and VectorDouble,
	// but make sure Python does not delete stuff.
	m.def("get_vector_int_ptr", &vector_int_to_cptr,
	      pybind11::return_value_policy::reference,
	      "Exposes the pointer to the raw data in a VectorInt");
	m.def("get_vector_double_ptr", &vector_double_to_cptr,
	      pybind11::return_value_policy::reference,
	      "Exposes the pointer to the raw data in a VectorInt");


	// Some options to add data fields:
	m.def("add_data_field", &lt_block_data_add_data_field,
	      "Adds a data field to given block_data.");
	m.def("add_special_field", &lt_block_data_add_special_field,
	      "Adds a data field as special field to given block_data.");

	m.def("set_meta", &lt_block_data_set_meta,
	      "Sets the meta-data of the block_data.");
	m.def("set_domain", &lt_block_data_set_domain,
	      "Sets the domain of the block_data.");

	m.def("get_domain_xlo", &lt_block_data_get_domain_xlo_vec,
	      "Returns the xlo of the domain.");
	m.def("get_domain_xhi", &lt_block_data_get_domain_xhi_vec,
	      "Returns the xhi of the domain.");

	m.def("get_domain_periodic", &lt_block_data_get_domain_periodic,
	      "Returns the periodic bits.");

	// And to mutate them:
	m.def("swap_fields", &lt_block_data_swap_fields);
	m.def("remove_field", &lt_block_data_remove_field);


	// Some debug stuff:
	m.def("print_stats", &lt_block_data_print_stats,
	      "Prints info about the C++-side of block_data");

	// To filter:
	m.def("filter_block_data", &lt_block_data_filter,
	      "Filters block_data based on indices." );

	// To create handles for new blocks:
	m.def( "new_block_data", &lt_new_block_data_handle,
	       "Creates a new block_data_handle." );
	m.def("delete_block_data", &lt_delete_block_data_handle,
	      "To delete a dynamically-created block_data_handle." );


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
		.value("IZ", lammps_tools::block_data::IZ)
		.value("ORIENT_X", lammps_tools::block_data::ORIENT_X)
		.value("ORIENT_Y", lammps_tools::block_data::ORIENT_Y)
		.value("ORIENT_Z", lammps_tools::block_data::ORIENT_Z)
		.value("ORIENT_W", lammps_tools::block_data::ORIENT_W);


	return m.ptr();
}
