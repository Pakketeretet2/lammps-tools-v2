#include "cluster_finder.hpp"
#include "id_map.hpp"
#include "util.hpp"

namespace lammps_tools {

void add_conns_to_network( neigh_list &conns,
                           std::list<std::list<int> > &networks )
{
	int max_mol = conns.size();
	std::vector<bool> mol_out( max_mol, false );

	for( int i = 1; i < mol_out.size(); ++i ){
		if( !mol_out[i] ){
			std::list<int> network;
			mol_out[i] = true;
			network.push_back( i );

			for( auto it = network.begin();
			     it != network.end(); ++it ){
				int j = *it;
				for( int k : conns[j] ){
					if( !mol_out[k] ){
						network.push_back(k);
						mol_out[k] = true;
					}
				}
			}
			networks.push_back( network );
		}
	}

}


void find_molecular_networks ( const block_data &b, const neigh_list &neighs,
                               neigh_list &conns,
                               std::list<std::list<int> > &networks )
{
	const std::vector<int> &mol = data_as<int>(
		b.get_special_field( block_data::MOL ) );
	const std::vector<int> &id = data_as<int>(
		b.get_special_field( block_data::ID ) );
	id_map im( id );
	int max_mol = *std::max_element( mol.begin(), mol.end() );
	conns.resize( max_mol + 1 );

	for( int i = 0; i < neighs.size(); ++i ){
		int mol_i = mol[i];
		for( int j : neighs[i] ){
			int mol_j = mol[j];

			if( mol_i == mol_j ) continue;

			if( util::contains( conns[mol_i], mol_j ) ||
			    util::contains( conns[mol_j], mol_i ) ) continue;

			conns[mol_i].push_back(mol_j);
			conns[mol_j].push_back(mol_i);
		}
	}

	// Reduce mol connections into networks:
	add_conns_to_network( conns, networks );
}


} // namespace lammps_tools
