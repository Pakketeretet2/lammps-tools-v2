#include <pybind11/pybind11.h>

#include "lt_block_writers.h"

PYBIND11_PLUGIN(block_writers_) {
	pybind11::module m("block_writers_",
	                   "Exposes functions to write block data" );

	m.def( "to_lammps_data", &lt_block_writers_lammps_data,
	       "Writes given block data to a LAMMPS data file." );
	m.def( "to_lammps_dump", &lt_block_writers_lammps_dump,
	       "Writes given block data to a LAMMPS dump file." );

	return m.ptr();
}
