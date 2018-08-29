#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "make_vectors_opaque.hpp"

/*
  Wrappers for C++ types:
*/
void print_list( pybind11::list my_list ) {
	for (auto item : my_list)
		std::cout << item << " ";
}

// This is probably dog slow:
std::vector<double> double_vec_from_list( pybind11::list my_list ) {
	std::vector<double> v;
	for (auto item : my_list){
		std::stringstream ss;
		ss << item;
		double value;
		ss >> value;
		v.push_back(value);
	}
	return v;
}

// This is probably dog slow:
std::vector<int> int_vec_from_list( pybind11::list my_list ) {
	std::vector<int> v;
	for (auto item : my_list){
		std::stringstream ss;
		ss << item;
		int value;
		ss >> value;
		v.push_back(value);
	}
	return v;
}





/*
std::vector<double> double_vec_from_np_arr( pybind11::array_t<double> arr )
{
	std::vector<double> v;
	for( auto x : arr ){
		std::stringstream ss;
		ss << x;
		double value;
		ss >> value;
		v.push_back( value );
	}
	return v;
}

std::vector<int> int_vec_from_np_arr( pybind11::array_t<int> arr )
{
	std::vector<int> v;
	for( auto x : arr ){
		std::stringstream ss;
		ss << x;
		int value;
		ss >> value;
		v.push_back( value );
	}
	return v;
}
*/


PYBIND11_PLUGIN(conversions_) {
	pybind11::module m("conversions_", "Converts Python --> C++");

	m.def("print_list", &print_list,
	      "Prints list from C++");

	m.def("list_to_double_vector", &double_vec_from_list,
	      "Converts list to double vector");
	m.def("list_to_int_vector", &int_vec_from_list,
	      "Converts list to int vector");

	/*
	m.def("numpy_array_to_double_vector", &double_vec_from_np_arr,
	      "Converts numpy array to double vector");
	m.def("numpy_array_to_int_vector", &int_vec_from_np_arr,
	      "Converts numpy array to int vector");
	*/
	return m.ptr();
}
