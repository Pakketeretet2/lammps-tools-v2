#ifndef NEIGHBORIZE_HPP
#define NEIGHBORIZE_HPP

/**
   \file neighborize.hpp

*/

#include "block_data.hpp"
#include "id_map.hpp"

#include <list>

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




/// Functor type for determining whether or not two atoms are neighbours.
struct are_neighbours
{
	virtual bool operator()( const block_data &b,
	                         int i, int j ) const = 0;
	virtual ~are_neighbours(){}
};


/// Standard distance criterion:
struct dist_criterion : public are_neighbours
{
	explicit dist_criterion( double rc ) : rc(rc), rc2(rc*rc) {}
	virtual bool operator()( const block_data &b, int i, int j ) const;
	virtual ~dist_criterion(){}
	double rc, rc2;
};


/**
   A shared interface for all neighborizers
*/
class neighborizer {
public:
	/// Sets how molecule status of atoms is taken into account.
	enum topological_policy {
		IGNORE = 0,  ///< Ignore molecule info
		INCLUDE = 1, ///< Always consider neighbours if in same mol
		EXCLUDE = -1 ///< Always consider not neighbours if in same mol.
	};

	/**
	   A generic neighbour list builder interface.

	   \note This expects that the special fields for both blocks are set.

	   \note This builds a neighbour list for all particles in
	         s1 w.r.t all particles in s2. That is, it finds for
	         all particles in s1 the neighbours in s2.

	   \param b1    block_data to neighborize.
	   \param s1    A list of indices to neighborise
	   \param s2    A list of indices to neighborise

	   \param dims  dimensionality of the system.
	*/

	neighborizer( const block_data &b, const std::list<int> &s1,
	              const std::list<int> &s2, int dims );

	virtual ~neighborizer(){}

	double build_list( neigh_list &neighs,
	                   const are_neighbours &criterion );


	int dims;
	int periodic;

	const block_data &b;

	double xlo[3];
	double xhi[3];

	int n_atoms;

	bool quiet;

	int mol_policy;
	int bond_policy;

	const std::list<int> &s1, &s2;
private:

	int append_particles_in_mol( neigh_list &neighs );
	int append_bonded_particles( neigh_list &neighs );

	int remove_particles_in_mol( neigh_list &neighs );

	virtual int build( neigh_list &neighs,
	                   const are_neighbours &criterion )
	{
		my_runtime_error(__FILE__, __LINE__, "virtual function called!");
		return 0;
	}

};

/**
   \brief Removes double entries in neigh list.

   \param neighs  neighbour list to clean.
*/
void remove_doubles( neigh_list &neighs );



/**
   \brief Calculates a neighbour list for the atoms in the block data.

   \param neighs[out]    The neighbour lists, indexed per particle ID
   \param b              block_data to neighborize.
   \param fields         Vector containing the data field names of id, type,
                         x, y, z, respectively.
   \param itype          Type of particle i to consider
   \param jtype          Type of particle j to consider
   \param method         Distance method to use (see neighborize_methods)
   \param dom            simulation domain
   \param dims           dimensions of the system (2 or 3)
   \param rc             Consider atoms less than this apart as neighbour
   \param mol_policy     Specifies how to take molecule ID into account.
   \param bond_policy    Specifies how to take bond topology into account.
   \param filter         Add only particles with filter(particle_ID) == true
   \param neigh_est      A guess for the number of neighbours.

   \returns The average number of neighbours per particle.
*/
double make_list_dist( neigh_list &neighs,
                       const block_data &b,
                       int itype, int jtype, int method, int dims,
                       double rc,
                       int mol_policy = neighborizer::IGNORE,
                       int bond_policy = neighborizer::IGNORE,
                       bool quiet = true );


/**
   \brief Calculates a neighbour list for the atoms in the block data.

   \param b              block_data to neighborize.
   \param fields         Vector containing the data field names of id, type,
                         x, y, z, respectively.
   \param itype          Type of particle i to consider
   \param jtype          Type of particle j to consider
   \param method         Distance method to use (see neighborize_methods)
   \param dom            simulation domain
   \param dims           dimensions of the system (2 or 3)
   \param rc             Consider atoms less than this apart as neighbour
   \param mol_policy     Specifies how to take molecule ID into account.
   \param bond_policy    Specifies how to take bond topology into account.
   \param filter         Add only particles with filter(particle_ID) == true
   \param neigh_est      A guess for the number of neighbours.

   \returns The average number of neighbours per particle.
*/
neigh_list nearest_neighs( const block_data &b,
                           int itype, int jtype, int method, int dims,
                           double rc,
                           int mol_policy = neighborizer::IGNORE,
                           int bond_policy = neighborizer::IGNORE,
                           bool quiet = true );


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
