#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

/*
  Wrappers for C++ types:
*/
std::vector<double> double_vec_from_pylist( pybind11::list list )
{
	std::vector<double> v;
	for( double x : list ){
		v.push_back(x);
	}
	return v;
}

std::vector<int> int_vec_from_pylist( pybind11::list list )
{
	std::vector<int> v;
	for( int x : list ){
		v.push_back(x);
	}
	return v;
}

std::vector<double> double_vec_from_np_arr( pybind11::array_t<double> arr )
{
	std::vector<double> v;
	for( double x : arr ){
		v.push_back(x);
	}
	return v;
}

std::vector<int> int_vec_from_np_arr( pybind11::array_t<int> arr )
{
	std::vector<int> v;
	for( int x : list ){
		v.push_back(x);
	}
	return v;
}

PYBIND11_PLUGIN(conversions_) {
	pybind11::module m("conversions_", "Converts Python --> C++");

	m.def("list_to_double_vector", &double_vec_from_pylist,
	      "Converts list to double vector");
	m.def("list_to_int_vector", &int_vec_from_pylist,
	      "Converts list to int vector");

	m.def("numpy_array_to_double_vector", &double_vec_from_np_arr,
	      "Converts numpy array to double vector");
	m.def("numpy_array_to_int_vector", &int_vec_from_np_arr,
	      "Converts numpy array to int vector");

	return m.ptr();
}
