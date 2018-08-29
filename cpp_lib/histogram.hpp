#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

/**
   \file histogram.hpp

   Contains routines for making histograms/distributions.
*/

#include <vector>

namespace lammps_tools {

namespace histogram {


/**
   \brief makes a histogram of given data.

   \param data    The data to histogram
   \param y0      Lower limit of values to consider
   \param y1      Upper limit of values to consider
   \param N_bins  Number of bins to use. Implicitly defines resolution

   \returns a (normalized) histogram of the data.
*/
template <typename T>
std::vector<double> make_histogram( const std::vector<T> &data,
                                    const T& y0, const T& y1, int N_bins );


/**
   \brief Specializes histogram for doubles.
*/
std::vector<double> make_histogram_double( const std::vector<double> &data,
                                           double y0, double y1, int N_bins );


/**
   \brief Specializes histogram for ints
*/
std::vector<double> make_histogram_int( const std::vector<int> &data,
                                        int y0, int y1, int N_bins );




std::vector<double> make_radial_grid( const std::vector<double> &xgrid,
                                      const std::vector<double> &ygrid,
                                      int Nbins );

template <typename T>
std::vector<T> radial_avg( const std::vector<double> &xgrid,
                           const std::vector<double> &ygrid,
                           const std::vector<T> &y, int Nbins );


std::vector<double> make_radial2_grid( const std::vector<double> &xgrid,
                                       const std::vector<double> &ygrid,
                                       int Nbins );

template <typename T>
std::vector<T> radial2_avg( const std::vector<double> &xgrid,
                            const std::vector<double> &ygrid,
                            const std::vector<T> &y, int Nbins );


std::vector<double> radial_avg( const std::vector<double> &xgrid,
                                const std::vector<double> &ygrid,
                                const std::vector<int> &y, int Nbins );


std::vector<double> radial_avg( const std::vector<double> &xgrid,
                                const std::vector<double> &ygrid,
                                const std::vector<int> &y, int Nbins );




} // namespace histogram

} // namespace lammps_tools



#endif // HISTOGRAM_HPP
