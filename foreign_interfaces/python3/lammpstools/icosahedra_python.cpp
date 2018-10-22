#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/icosahedra.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)

PYBIND11_PLUGIN(icosahedra_) {
	using namespace lammps_tools;

	pybind11::module m("icosahedra", "Exposes icosahedra analysis.");

	m.def( "dihedral_angles", &icosahedra::get_dihedral_angles,
	       "Get the dihedral angles between triangles." );

	m.def( "get_mol_connections", &icosahedra::make_mol_connections,
	       "Generates the neighbor list between molecules." );

	m.def( "patch_types_to_ids", &icosahedra::patch_types_to_ids,
	       "Converts patch types to an id list. ");

	m.def( "get_scattering_intensity", &icosahedra::get_scattering_intensity,
	       "Calculates an effective, relative scattering intensity." );


	return m.ptr();
}
