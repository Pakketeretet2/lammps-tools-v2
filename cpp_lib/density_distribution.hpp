#ifndef DENSITY_DISTRIBUTION_HPP
#define DENSITY_DISTRIBUTION_HPP

/**
   \file density_distribution.hpp

   Contains functions that analyze the distribution of particles in the box.
*/

#include "block_data.hpp"

#include <cmath>

namespace lammps_tools {

namespace density {


/**
   \brief "boxes" atoms in a grid

   \param b     Block data to box
   \param Nx    Number of bins to make.
   \param dx    Will contain the grid spacing
   \param dims  Dimension of systems

   \returns a vector that maps atom index to box index.
*/
std::vector<int> box_atoms( const lammps_tools::block_data &b,
                            int Nx, double &dx, int dims );


/**
   \brief constructs a density distribution on a grid.

   \param b     Block data to box
   \param Nx    Number of bins to make
   \param itype Atom type to count
   \param dims  Dimension of systems
*/
std::vector<double> density_distribution( const lammps_tools::block_data &b,
                                          int Nx, int itype, int dims );



/**
   \brief calculates the number of atoms in each box

   \param b       Block data to analyze
   \param boxes   Vector that specifies the box per atom.
   \param N_boxes The number of boxes there are.
   \param itype   The atom type to count
   \param quiet   If false, print some output on progress.

   \returns a vector of counts per box.
*/
std::vector<int> box_count( const lammps_tools::block_data &b,
                            const std::vector<int> &boxes,
                            int N_boxes, int itype, bool quiet );


/**
   \brief Calculates the distance between two boxes.

   \param bin_i    The index of box i
   \param bin_j    The index of box j
   \param dx       The box width
   \param Nx       The number of boxes in each dim.
   \param dim      The number of dimensions.
*/
double box_distance( int bin_i, int bin_j, double dx,
                     int Nx, int dim );


} // namespace density

} // namespace lammps_tools




#endif // DENSITY_DISTRIBUTION_HPP
