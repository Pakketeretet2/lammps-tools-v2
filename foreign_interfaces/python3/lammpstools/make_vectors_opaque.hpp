#ifndef MAKE_VECTORS_OPAQUE_HPP
#define MAKE_VECTORS_OPAQUE_HPP

#include <complex>
#include <vector>

PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)
PYBIND11_MAKE_OPAQUE(std::vector<std::complex<int> >)
PYBIND11_MAKE_OPAQUE(std::vector<std::complex<double> >)


#endif // MAKE_VECTORS_OPAQUE_HPP
