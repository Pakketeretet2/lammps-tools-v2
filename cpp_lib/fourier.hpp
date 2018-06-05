#ifndef FOURIER_HPP
#define FOURIER_HPP

#include <complex>
#include <vector>


namespace lammps_tools {

namespace fourier {

typedef std::complex<double> cx_double;


std::vector<double> fft_abs( const std::vector<cx_double > &fft );

std::vector<double> fft_real( const std::vector<cx_double > &fft );

std::vector<double> fft_imag( const std::vector<cx_double > &fft );



std::vector<cx_double > fft( int Nx, int Ny, int Nz,
                             const std::vector<double> &data );

std::vector<cx_double > fft_shift( int  Nx, int Ny, int Nz,
                                   const std::vector<cx_double > &fft );

std::vector<cx_double > ifft_shift( int  Nx, int Ny, int Nz,
                                    const std::vector<cx_double > &fft );


std::vector<cx_double > fft_double( int Nx, int Ny, int Nz,
                                    const std::vector<double> &data );

std::vector<cx_double > fft_int( int Nx, int Ny, int Nz,
                                 const std::vector<int> &data );






} // namespace fourier

} // lammps_tools



#endif // FOURIER_HPP
