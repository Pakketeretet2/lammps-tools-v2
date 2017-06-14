#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lt_block_data.h"
#include "lt_transformations.h"

PYBIND11_PLUGIN(transformations_) {
	pybind11::module m("transformations_",
	                   "Exposes transformations through pybind11" );

	m.def("rotate_all", &lt_transformations_rotate_all_pb,
	      "Rotates all particles over an axis at given origin.");
	m.def("shift_all", &lt_transformations_shift_all_pb,
	      "Shifts all particles over given distance.");
	/*m.def("center_box_on", &lt_transformations_center_box_on_pb,
	      "Centers the entire box on given point.");
*/
	m.def("shift", &lt_transformations_shift_pb,
	      "Shifts given particles over given distance");
	m.def("rotate", &lt_transformations_rotate_pb,
	      "Rotates given particles over an axis at given origin.");


	return m.ptr();
}
