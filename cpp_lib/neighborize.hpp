#ifndef NEIGHBORIZE_HPP
#define NEIGHBORIZE_HPP

/**
   \file neighborize.hpp

*/

#include "block_data.hpp"

namespace lammps_tools {

namespace neighborize {


/// Typedef for a neigh list:
typedef std::vector<std::vector<int> > neigh_list;


/**
  Enum of the neighborization methods.
*/
enum neighborize_methods {
	DIST_NSQ    = 0, ///< A distance criterion, using a method slow for big data sets
	DIST_BIN    = 1, ///< A distance criterion, using a method fast for big data sets
        DELAUNAY    = 2, ///< Delaunay triangulation in 2/3D, using the CGAL lib
	CONVEX_HULL = 3  ///< Convex hull for particles on ellipsoids/spheres, using CGAL lib
};



/// A typedef for a function that can be used to filter particles.
typedef bool( *particle_filter )( const block_data &, int );


/// The default filter is to let everything pass:
inline bool pass_all( const block_data &, int ) { return true; }


/// Functor type for determining whether or not two atoms are neighbours.
struct are_neighbours
{
	virtual bool operator()( const block_data &b, int i, int j ) const = 0;
};


/// Standard distance criterion:
struct dist_criterion : public are_neighbours
{
	dist_criterion( double rc ) : rc(rc), rc2(rc*rc) {}
	virtual bool operator()( const block_data &b, int i, int j ) const;
	double rc, rc2;
};



/**
   \brief Calculates a neighbour list for the atoms in the block data.
   
   \param neighs[out]    The neighbour lists, indexed per particle_index
   \param b              block_data to neighborize.
   \param fields         Vector containing the data field names of id, type,
                         x, y, z, respectively.
   \param itype          Type of particle i to consider
   \param jtype          Type of particle j to consider
   \param method         Distance method to use (see neighborize_methods)
   \param dom            simulation domain
   \param dims           dimensions of the system (2 or 3)
   \param rc             Consider atoms less than this apart as neighbour
   \param include_mol    Count every atom in molecule as neighbour
   \param include_bonds  Count every atom in molecule as neighbour
   \param filter         Add only particles with filter(particle_index) == true
   \param neigh_est      A guess for the number of neighbours.

   \returns The average number of neighbours per particle.
*/
double make_list_dist( neigh_list &neighs,
                       const block_data &b,
                       const std::vector<std::string> &fields,
                       int itype, int jtype, int method, int dims,
                       double rc,
                       bool include_mol = false,
                       bool include_bonds = false,
                       particle_filter filt = pass_all,
                       int neigh_est = 24 );


/**
   \brief Verifies neigh list.

   It checks
      1. if i is neigh of j, j is neigh of i.
      2. no double counting.

   Raises a logic_error if the check does not check out.
*/
void verify_neigh_list( const neigh_list &neighs );

} // namespace neighborize

} // namespace lammps_tools


#endif // NEIGHBORIZE_HPP
