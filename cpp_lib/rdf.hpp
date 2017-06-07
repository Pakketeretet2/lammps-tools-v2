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

void compute_rdf( const block_data &b, int Nbins, double r0, double r1,
                  int dims, int itype, int jtype,
                  std::vector<double> &rdf, std::vector<double> &coord );

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
