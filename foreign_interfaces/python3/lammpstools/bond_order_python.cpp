#include "pybind11/pybind11.h"
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "lt_block_data.h"

#include "../../../cpp_lib/bond_order.hpp"

#include "make_vectors_opaque.hpp"

// Converts a numpy array to a vector so that it can be passed to bond_order
namespace {

using namespace lammps_tools;
//using c_style = pybind11::array::c_style;

double compute_psi_n( const block_data &b,
                      const neighborize::neigh_list &neighs,
                      int n,
                      pybind11::array_t<double> &aaxis,
                      pybind11::array_t<double> &psi_n_real,
                      pybind11::array_t<double> &psi_n_imag )
{
	std::vector<double> aangles( b.N );
	auto buff_axis = aaxis.request();
	double *axis_arr = static_cast<double*>( buff_axis.ptr );
	point axis( axis_arr[0], axis_arr[1], axis_arr[2] );

	std::vector<double> psi_n_real_v( b.N );
	std::vector<double> psi_n_imag_v( b.N );

	auto buff_real = psi_n_real.request();
	auto buff_imag = psi_n_imag.request();
	double *psi_n_real_arr = static_cast<double*>( buff_real.ptr );
	double *psi_n_imag_arr = static_cast<double*>( buff_imag.ptr );

	double val = order_parameters::compute_psi_n( b, neighs, n, axis,
	                                              psi_n_real_v,
	                                              psi_n_imag_v );
	for( std::size_t i = 0; i < psi_n_real_v.size(); ++i ){
		psi_n_real_arr[i] = psi_n_real_v[i];
		psi_n_imag_arr[i] = psi_n_imag_v[i];
	}

	return val;
}

}


PYBIND11_PLUGIN(bond_order_) {
	pybind11::module m("block_data_", "Exposes bond order functions");
	using namespace lammps_tools;

	m.def( "compute_psi_n", &compute_psi_n,
	       "Calculates the bond order parameter." );


	return m.ptr();
}
