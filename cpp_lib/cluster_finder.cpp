#include "cluster_finder.hpp"
#include "id_map.hpp"
#include "util.hpp"

namespace lammps_tools {

namespace neighborize {

void add_conns_to_network( const neigh_list &conns, neigh_list &networks )
{
	int max_mol = conns.size();
	std::vector<bool> mol_out( max_mol, false );
	for( std::size_t i = 1; i < mol_out.size(); ++i ){
		if( !mol_out[i] ){
			std::vector<int> network;
			mol_out[i] = true;
			network.push_back( i );

			for( std::size_t it = 0; it < network.size(); ++it ){
				int j = network[it];
				const std::vector<int> &j_neighs = conns[j];
				for( int o_mol : j_neighs ){
					if( mol_out[o_mol] ) continue;

					mol_out[o_mol] = true;
					network.push_back(o_mol);
				}
			}

			networks.push_back( network );
		}
	}

}



neigh_list get_molecular_connections( const block_data &b,
                                      const neigh_list &neighs )
{
	const std::vector<int> &mol = data_as<int>(
		b.get_special_field( block_data::MOL ) );

	int max_mol = *std::max_element( mol.begin(), mol.end() );
	neigh_list conns( max_mol + 1 );
	const std::vector<int> &type = data_as<int>(
		b.get_special_field( block_data::TYPE ) );

	for( std::size_t i = 0; i < neighs.size(); ++i ){
		int mol_i = mol[i];
		for( int j : neighs[mol_i] ){
			int mol_j = mol[j];
			if( mol_i == mol_j ) continue;

			if( util::contains( conns[mol_i], mol_j ) ||
			    util::contains( conns[mol_j], mol_i ) ) continue;

			conns[mol_i].push_back(mol_j);
			conns[mol_j].push_back(mol_i);
		}
	}

	return conns;
}



} // namespace neighborize

} // namespace lammps_tools
