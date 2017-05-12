#include <catch.hpp>

#include <iostream>
#include <iomanip>

#include "block_data.hpp"
#include "data_field.hpp"
#include "util.hpp"

TEST_CASE( "Sorting works", "[util_sort]" ) {
	using namespace lammps_tools;

	auto double_comp = []( double d1, double d2 ){ return d1 < d2; };
	std::vector<double> a = { 2.0, 1.0, 5.2, 1.2, 14.0 -23.0 };
	std::vector<std::size_t> p = util::sort_permutation( a, double_comp );

	std::cerr << "Before:\n";
	for( double x : a ){
		std::cerr << "  " << x << "\n";
	}
	std::cerr << "Permutation:\n";
	for( std::size_t x : p ){
		std::cerr << "  " << x << "\n";
	}
	std::vector<double>      target_vec  = { -9, 1, 1.2, 2, 5.2 };
	std::vector<std::size_t> target_perm = { 4,1,3,0,2 };

	std::cerr << "Sort should be:\n";
	for( std::size_t i = 0; i < p.size(); ++i ){
		std::cerr << "  " << a[p[i]] << "\n";
		REQUIRE( p[i] == target_perm[i] );
	}
	util::apply_permutation( a, p );
	std::cerr << "Sort is:\n";
	for( std::size_t i = 0; i < a.size(); ++i ){
		std::cerr << "  " << a[i] << "\n";
		REQUIRE( a[i] == target_vec[i] );
	}
}


TEST_CASE( "Sorting block works", "[block_data_sort]" ) {
	using namespace lammps_tools;
	data_field_double d1( "d1", 5 );
	data_field_double d2( "d2", 5 );
	data_field_int    d3( "d3", 5 );
	data_field_double d4( "d4", 5 );

	d1[0] = 3.0;
	d1[1] = 2.2;
	d1[2] = 12.1;
	d1[3] = -3.2;
	d1[4] = -4.1;

	d2[0] = 2.0;
	d2[1] = -2.2;
	d2[2] = -12.3;
	d2[3] =  4.2;
	d2[4] =  4.5;

	d3[0] = 1;
	d3[1] = 2;
	d3[2] = 3;
	d3[3] = 4;
	d3[4] = 5;

	d4[0] = 1.0;
	d4[1] = 0.5;
	d4[2] = 0.0;
	d4[3] = 0.5;
	d4[4] = 1.0;

	data_field_double *d1o = static_cast<data_field_double*>( copy(&d1) );
	data_field_double *d2o = static_cast<data_field_double*>( copy(&d2) );
	data_field_int    *d3o = static_cast<data_field_int*>( copy(&d3) );
	data_field_double *d4o = static_cast<data_field_double*>( copy(&d4) );


	SECTION( "Sort data fields directly." ){
		auto comp = []( double a, double b ) { return a < b; };
		std::vector<std::size_t> p =
			get_data_permutation<double>( &d1, comp );
		std::cerr << "Before sort:\n";
		for( int i = 0; i < 5; ++i ){
			std::cerr << "    " << std::setw(5) <<  d1[i]
			          << "\t" << std::setw(5) << d2[i]
			          << "\t" << std::setw(5) << d3[i]
			          << "\t" << std::setw(5) << d4[i] << "\n";
		}
		std::cerr << "Permutation:";
		for( std::size_t i : p ){
			std::cerr << " " << i;
		}
		std::cerr << "\nAfter sort:\n";
		sort_data_with_permutation<double>( &d1, p );
		sort_data_with_permutation<double>( &d2, p );
		sort_data_with_permutation<int>( &d3, p );
		sort_data_with_permutation<double>( &d4, p );


		for( int i = 0; i < 5; ++i ){
			std::cerr << "    " << std::setw(5) <<  d1[i]
			          << "\t" << std::setw(5) << d2[i]
			          << "\t" << std::setw(5) << d3[i]
			          << "\t" << std::setw(5) << d4[i] << "\n";
			REQUIRE( d1[i] == (*d1o)[p[i]] );
			REQUIRE( d2[i] == (*d2o)[p[i]] );
			REQUIRE( d3[i] == (*d3o)[p[i]] );
			REQUIRE( d4[i] == (*d4o)[p[i]] );
		}
		for( int i = 0; i < 4; ++i ){
			REQUIRE( d1[i+1] >= d1[i] );
		}
	}

	SECTION( "Sort data fields in block_data." ){
		block_data b;
		using v_d = std::vector<double>;
		using v_i = std::vector<int>;
		b.set_natoms( 5 );

		std::vector<std::string> fields = { "d1", "d2", "d3", "d4" };
		b.add_field( d1 );
		b.add_field( d2 );
		b.add_field( d3 );
		b.add_field( d4 );

		std::cerr << "Before sort:\n";
		for( int i = 0; i < 5; ++i ){
			std::cerr << "  ";
			for( const std::string &n : fields ){
				const data_field *df = b.get_data(n);
				int type = df->type();
				std::cerr << "\t" << std::setw(5);
				if( type == data_field::INT ){
					const v_i& v = data_as<int>( df );
					std::cerr  << v[i];
				}else if( type == data_field::DOUBLE ){
					const v_d& v = data_as<double>( df );
					std::cerr << v[i];
				}
			}
			std::cerr << "\n";
		}

		auto comp = []( double a, double b ) { return a < b; };
		b.sort_along( "d1", comp );

		std::cerr << "After sort:\n";
		for( int i = 0; i < 5; ++i ){
			std::cerr << "  ";
			for( const std::string &n : fields ){
				const data_field *df = b.get_data(n);
				int type = df->type();
				std::cerr << "\t" << std::setw(5);
				if( type == data_field::INT ){
					const v_i& v = data_as<int>( df );
					std::cerr  << v[i];
				}else if( type == data_field::DOUBLE ){
					const v_d& v = data_as<double>( df );
					std::cerr << v[i];
				}
			}
			std::cerr << "\n";
		}

	}

	delete d1o;
	delete d2o;
	delete d3o;
	delete d4o;

}
