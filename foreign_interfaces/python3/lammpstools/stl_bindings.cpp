#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);


PYBIND11_PLUGIN (stl_binds) {
	pybind11::module m("stl_binds", "Dummy module that binds STL containers opaquely.");

	pybind11::bind_vector<std::vector<int>>(m, "VectorInt");
	pybind11::bind_vector<std::vector<double>>(m, "VectorDouble");


}
