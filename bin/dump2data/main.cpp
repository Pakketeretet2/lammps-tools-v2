#include <iostream>
#include <memory>
#include <vector>

#include "data_field.hpp"
#include "block_data_access.hpp"
#include "writers_lammps.hpp"
#include "dump_reader_lammps.hpp"
#include "util.hpp"

void print_usage()
{
	std::cerr << "Use this like dump2data <dump file> -c <column headers> "
	          << " -o <output>,\n"
	          << "  where <headers> is a quoted string containing the\n"
	          << "  names of the headers in the dump file and <output> is\n"
	          << "  an optional output dump file name.\n";
}


int main( int argc, char **argv )
{
	using namespace lammps_tools;

	std::string dump;
	std::string headers = "";
	std::string out_file = "-";
	bool to_local = false;
	bool silent = false;
	bool ignore_first = false;
	bool recreate_img_flags = true;
	int hack_mol_stride = 0;

	if (argc < 2) {
		std::cerr << "Pass a dump file!\n";
		return -2;
	}

	dump = argv[1];
	if (dump == "-h" || dump == "--help") {
		// not actually a dump but user wants help:
		print_usage();
		return 1;
	}

	int i = 2;
	while( i < argc ){
		std::string a( argv[i] );
		if( a[0] == '-' ){
			if( a == "-c" || a == "--column-headers" ){
				headers = argv[i+1];
				i += 2;
			}else if( a == "-o" || a == "--output" ){
				out_file = argv[i+1];
				i += 2;
			}else if( a == "-l" || a == "--local" ){
				to_local = true;
				i += 1;
			}else if( a == "-h" || a == "--help" ){
				print_usage();
				return 1;
			}else if( a == "-i" || a == "--ignore-first" ){
				ignore_first = true;
				i += 1;
			}else if( a == "-s" || a == "--silent" ){
				silent = true;
				i += 1;
			}else if (a == "-r" || a == "--no-recreate-img-flags") {
				recreate_img_flags = false;
				i += 1;
			} else if( a == "-m" || a == "--hack-mol-stride" ) {
				hack_mol_stride = std::stoi(argv[i+1]);
				i += 2;
			}
		} else {
			std::cerr << "Unrecognized trailing argument \""
			          << a << "!\n";
			return -1;
		}
	}

	if (hack_mol_stride < 0) {
		std::cerr << "hack_mol_stride cannot be negative!\n";
		return -7;
	}

	int fformat = FILE_FORMAT_PLAIN;
	if (util::ends_with(dump, ".bin")) {
		fformat = FILE_FORMAT_BIN;
	} else if (util::ends_with(dump, ".gz")) {
		fformat = FILE_FORMAT_GZIP;
	}


	if (headers.empty() && (fformat == FILE_FORMAT_BIN)) {
		std::cerr << "You need to provide column headesr for binary dump files!\n";
		return -6;
	}

	std::cerr << "Headers are " << headers << ".\n";

	std::ofstream *out_fstream = nullptr;
	std::ostream *out = nullptr;
	int out_fformat = FILE_FORMAT_PLAIN;

	if( out_file == "-" ){
		out = &std::cout;
	}else{
		if( util::ends_with( out_file, ".bin") ){
			std::cerr << "Data files can only be binary!\n";
			return -3;
		}else if( util::ends_with( out_file, ".gz" ) ){
			std::cerr << "Writing to gzip is not supported by "
			          << "dump-cat. Instead, stream text output "
			          << "through gzip!\n";
			return -4;
		}else{
			std::cerr << "Assuming file " << out_file
			          << " is plain text dump file.\n";
			out_fformat = FILE_FORMAT_PLAIN;
			out_fstream = new std::ofstream( out_file );
			out = out_fstream;
		}
	}

	int n_writes = 0;
	auto status_print = [&n_writes, silent]	{
		if( !silent && n_writes > 0  && n_writes % 50 == 0 ){
			std::cerr << "At block " << n_writes << "...\n";
		}
		n_writes++; };

	block_data b;
	std::unique_ptr<readers::dump_reader_lammps> d(
	                                               readers::make_dump_reader_lammps( dump, fformat ) );
	std::vector<std::string> headers_vec = util::split(headers);
	d->set_column_headers( headers_vec );

	for (std::string h : headers_vec) {
		if (h == "id" || h == "ID") {
			d->set_column_header_as_special(h, block_data::ID);
		} else if (h == "type" || h == "TYPE") {
			d->set_column_header_as_special(h, block_data::TYPE);
		} else if (h == "x" || h == "X") {
			d->set_column_header_as_special(h, block_data::X);
		} else if (h == "y" || h == "Y") {
			d->set_column_header_as_special(h, block_data::Y);
		} else if (h == "z" || h == "Z") {
			d->set_column_header_as_special(h, block_data::Z);
		} else if (h == "mol" || h == "MOL") {
			d->set_column_header_as_special(h, block_data::MOL);
		}

	}

	while( d->next_block(b) == 0 ){
		//writers::block_to_lammps_dump( *out, b,
		//                               out_fformat,
		//                               to_local );
		status_print();
	}

	if (hack_mol_stride) {
		// Create a new mol column and use the stride to deduce mol.
		std::vector<int> mol(b.N);
		auto ids = get_id(b);
		for (long i = 0; i < b.N; ++i) {
			int idm = ids[i] - 1;
			mol[i] = 1 + idm/hack_mol_stride;
		}

		data_field_int mol_field( "mol", mol );
		b.add_field(mol_field, block_data::MOL);
	}

	if (recreate_img_flags) {
		std::vector<int> imx(b.N);
		std::vector<int> imy(b.N);
		std::vector<int> imz(b.N);

		b.dom.reconstruct_image_flags(b, imx, imy, imz);
		data_field_int img_x("ix", imx);
		data_field_int img_y("iy", imy);
		data_field_int img_z("iz", imz);

		b.add_field(img_x, block_data::IX);
		b.add_field(img_y, block_data::IY);
		b.add_field(img_z, block_data::IZ);

	}

	writers::block_to_lammps_data(*out, b);
	if( out_fstream ) delete out_fstream;
}
