#ifndef CORRELATION_HPP
#define CORRELATION_HPP

/**
   \file correlation.hpp

   This file contains basic routines for correlating data.
*/
#include "block_data.hpp"
#include "data_field.hpp"


namespace lammps_tools {

namespace correlate {


/**
   \brief Correlate a value over distance.

   Specifically, if the given data is f, then this calculates
   <f(r)*f(r+R)> - <f(r)>*<f(r+R)>, where the <> indicates averaging
   over all atoms.

   \param b    The block data to calculate the correlation for
   \param data The data field to average
   \param x0   Lower point on distance over which to correlate
   \param x1   Higher point on distance over which to correlate
   \param dx   Bin size.
   \param dims Dimensions of the simulation box

   \returns the correlation of given data.
*/
std::vector<double> correlate_int( const lammps_tools::block_data &b,
                                   const std::vector<int> &data,
                                   double x0, double x1, double dx, int dims );


std::vector<double> correlate_double( const lammps_tools::block_data &b,
                                      const std::vector<double> &data,
                                      double x0, double x1, double dx, int dims );






} // correlate

} // namespace lammps_tools


#endif // CORRELATION_HPP
