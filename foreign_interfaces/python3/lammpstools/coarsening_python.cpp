#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../../cpp_lib/coarsening.hpp"

#include "make_vectors_opaque.hpp"


PYBIND11_PLUGIN(coarsening_) {
	using namespace lammps_tools;

	pybind11::module m("coarsening_", "Exposes coarsening functions.");

	//m.def( "atoms_to_box", &coarse::atoms_to_box,
	//       "Converts atom positions to boxes." );

	m.def( "box_to_psi", &coarse::box_to_psi,
	       "Converts box occupation to psi order parameter." );

	m.def( "map_to_atom_double", &coarse::map_to_atom_double,
	       "Converts per box data back to psi per atom." );

	m.def( "map_to_atom_int", &coarse::map_to_atom_int,
	       "Converts per box data back to psi per atom." );

	m.def( "map_to_atom_cx_double", &coarse::map_to_atom_cx_double,
	       "Converts per box data back to psi per atom." );

	m.def( "map_to_atom_cx_int", &coarse::map_to_atom_cx_int,
	       "Converts per box data back to psi per atom." );


	return m.ptr();
}
