#ifndef COARSENING_HPP
#define COARSENING_HPP


/**
   \file coarsening.hpp Contains routines to quantify coarsening of mixtures

*/


#include "block_data.hpp"
#include "id_map.hpp"

#include <complex>
#include <vector>

namespace lammps_tools {

namespace coarse {


typedef std::complex<int> cx_int;
typedef std::complex<double> cx_double;



/**
   \brief calculates the psi value for each atom.

   \param b       Block data to analyze
   \param boxes   Vector that specifies the box per atom.
   \param N_boxes The number of boxes there are.
   \param seed    Random seed used to assign a psi for na==nb.
   \param quiet   If false, print some output on progress.

   \returns a vector of psis per box.
*/
std::vector<int> box_to_psi( const lammps_tools::block_data &b,
                             const std::vector<int> &boxes,
                             int N_boxes, int seed, bool quiet );



/**
   \brief Maps given per-box data back to all atoms.

   \param b     Block data to analyze
   \param boxes Specifies the box per atom
   \param data  The per-box info that is to be mapped onto atoms.

   \returns the data per atom.
*/
std::vector<double> map_to_atom_double( const lammps_tools::block_data &b,
                                        const std::vector<int> &boxes,
                                        const std::vector<double> &data );

/**
   \brief Maps given per-box data back to all atoms.

   \param b     Block data to analyze
   \param boxes Specifies the box per atom
   \param data  The per-box info that is to be mapped onto atoms.

   \returns the data per atom.
*/
std::vector<int> map_to_atom_int( const lammps_tools::block_data &b,
                                  const std::vector<int> &boxes,
                                  const std::vector<int> &data );


/**
   \brief Maps given per-box data back to all atoms.

   \param b     Block data to analyze
   \param boxes Specifies the box per atom
   \param data  The per-box info that is to be mapped onto atoms.

   \returns the data per atom.
*/
std::vector<cx_double> map_to_atom_cx_double( const lammps_tools::block_data &b,
                                              const std::vector<int> &boxes,
                                              const std::vector<cx_double> &data );



/**
   \brief Maps given per-box data back to all atoms.

   \param b     Block data to analyze
   \param boxes Specifies the box per atom
   \param data  The per-box info that is to be mapped onto atoms.

   \returns the data per atom.
*/
std::vector<cx_int> map_to_atom_cx_int( const lammps_tools::block_data &b,
                                        const std::vector<int> &boxes,
                                        const std::vector<cx_int> &data );





/**
   \brief Calculates the Fourier transform of given boxed data.

   \param Nx   Number of x positions of grid
   \param Ny   Number of y positions of grid
   \param Nz   Number of z positions of grid
   \param data Data to calculate the Fourier transform of.

   \returns The Fourier transform of given data.
*/
std::vector<double> fourier_transform( int Nx, int Ny, int Nz,
                                       const std::vector<double> &data );



} // namespace coarse

} // lammps_tools


#endif // COARSENING_HPP
