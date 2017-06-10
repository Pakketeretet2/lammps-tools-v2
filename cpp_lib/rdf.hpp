#ifndef RDF_HPP
#define RDF_HPP

/**
  \file rdf.hpp

  \brief Contains routines to calculate pair correlation functions and the like.
*/

#include <vector>

#include "block_data.hpp"
#include "neighborize.hpp"

namespace lammps_tools {

namespace neighborize {

/**
   \brief Calculate g(r) for itype with respect to jtype.

   \param[in]  b      The block_data to determine g(r) for.
   \param[in]  Nbins  Bin resolution
   \param[in]  r0     Lower cutoff for distance
   \param[in]  r1     Upper cutoff for distance
   \param[in]  dims   Dimensions of the system (used for normalisation)
   \param[in]  itype  First type
   \param[in]  jtype  Second type
   \param[out] rdf    Will contain the rdf
   \param[out] coord  Will contain the coordination (integral of rdf)

   \note This function constructs an internal neighbour list. If a
   neighbour list that ranges to at least r1 is available, consider
   calling compute_rdf_with_neighs for efficiency.
*/
void compute_rdf( const block_data &b, int Nbins, double r0, double r1,
                  int dims, int itype, int jtype,
                  std::vector<double> &rdf, std::vector<double> &coord );

/**
   \brief Calculate g(r) for itype with respect to jtype.

   \param[in]  b      The block_data to determine g(r) for.
   \param[in]  Nbins  Bin resolution
   \param[in]  r0     Lower cutoff for distance
   \param[in]  r1     Upper cutoff for distance
   \param[in]  dims   Dimensions of the system (used for normalisation)
   \param[in]  neighs Neighbour list to use for calculating the rdf
   \param[out] rdf    Will contain the rdf
   \param[out] coord  Will contain the coordination (integral of rdf)
*/
void compute_rdf_with_neighs( const block_data &b, int Nbins,
                              double r0, double r1, int dims,
                              const neigh_list &neighs,
                              std::vector<double> &rdf,
                              std::vector<double> &coord );

std::vector<double> rdf( const block_data &b, int Nbins, double r0, double r1,
                         int dims, int itype, int jtype );

std::vector<double> coord( const block_data &b, double r0, double r1, int dims,
                           const std::vector<double> &rrdf );






} // namespace neighborize

} // namespace lammps_tools


#endif // RDF_HPP
