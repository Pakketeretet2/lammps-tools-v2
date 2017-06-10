#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../../cpp_lib/block_data.hpp"
#include "../../../cpp_lib/scatter.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)

PYBIND11_PLUGIN(scatter_) {
	using namespace lammps_tools;
	pybind11::module m("scatter", "Functions to simulate scattering data.");

	m.def("rayleigh_gans_", &scatter::rayleigh_gans_,
	      "Calculates Rayleigh-Gans scattering data for specific ids.");

	return m.ptr();
}
