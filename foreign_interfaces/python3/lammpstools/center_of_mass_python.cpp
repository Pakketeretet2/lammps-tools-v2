#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../../cpp_lib/block_data.hpp"
#include "../../../cpp_lib/center_of_mass.hpp"


std::vector<double> center_of_mass_vec( const lammps_tools::block_data &b )
{
	lammps_tools::point p = lammps_tools::center_of_mass( b );
	std::vector<double> x(3);
	x[0] = p.x;
	x[1] = p.y;
	x[2] = p.z;
	return x;
}

std::vector<double> geometric_center_vec( const lammps_tools::block_data &b )
{
	lammps_tools::point p = lammps_tools::geometric_center( b );
	std::vector<double> x(3);
	x[0] = p.x;
	x[1] = p.y;
	x[2] = p.z;
	return x;
}



PYBIND11_PLUGIN(center_of_mass_) {

	pybind11::module m("center_of_mass_", "Calculates center of mass");

	m.def( "center_of_mass",   &center_of_mass_vec );
	m.def( "geometric_center", &geometric_center_vec );

	return m.ptr();
}
