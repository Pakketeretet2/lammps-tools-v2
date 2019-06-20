#include <iomanip>
#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <sstream>

#include "readers.hpp"
#include "dump_reader_lammps_bin.hpp"
#include "my_timer.hpp"
#include "util.hpp"
#include "writers.hpp"

using namespace lammps_tools;


// Simple program to extract frames from a dump file.
void print_usage()
{
	std::cerr << "Usage: extract-dump <dump-file> <options>\n"
	          << "Possible options:\n"
	          << "    -f/--file-format:  Specify dump file format (plain, bin, gzip)\n"
	          << "    -i/--info:         Print info about dump file\n"
	          << "    -n/--nframes:      A list of which frames to dump\n"
	          << "    -c/--columns:      Specifies the column headers.\n";
}




void print_dump_info(std::string dump_file, int file_format,
                     const std::vector<std::string> &headers)
{
	using namespace lammps_tools::readers;
	std::unique_ptr<dump_reader> dr(make_dump_reader(dump_file,
	                                                 file_format,
	                                                 DUMP_FORMAT_LAMMPS));
	dr->set_column_headers(headers);
	int n_frames = 0;
	block_data b;
	my_timer timer;
	double performance = 0;
	double elapsed = 0;
	while (dr->next_block(b) == 0) {
		++n_frames;
		if (n_frames % 50 == 0) {
			elapsed = timer.toc();
			performance = 1000*n_frames / elapsed;
			std::cerr << "  At frame " << n_frames << " ("
			          << performance << " blocks/s)...\n";
		}
	}
	elapsed = timer.toc();
	std::cerr << "Reading out all frames took " << elapsed/1000 << " s ("
	          << elapsed << " ms). Performance: " << performance
	          << " blocks/s.\n";
	std::cout << "Dump file:         " << dump_file << "\n";
	std::cout << "Number of frames:  " << n_frames  << "\n";
}



// Reads a vector from argv, starting at i.
// i is modified to point to one past the end of the vector.
template <typename T>
std::vector<T> read_vector(int argc, char **argv, int &i)
{
	std::vector<T> vec;
	std::string current = argv[i];
	std::cerr << "Reading vector... ";
	if (current != "{") {
		std::cerr << "Expected a brace-enclosed initializer list!\n";
		return vec;
	}
	std::cerr << "{";
	++i;
	while (i < argc) {
		current = argv[i++];
		if (current == "}") {
			std::cerr << " }\n";
			return vec;
		}

		std::stringstream ss;
		ss << current;
		T value;
		ss >> value;
		std::cerr << " " << value;
		vec.push_back(value);
	}
	// If you reached this you did not encounter the closing brace.
	// Should this be an error?
	std::cerr << "\n";
	std::cerr << "Did not find closing brace \"}\"!\n";
	return vec;
}



int main(int argc, char **argv)
{
	int i = 2;
	int file_format = FILE_FORMAT_PLAIN;
	std::string dump_file = argv[1];
	bool dump_info = false;
	std::vector<int> frames;
	std::vector<std::string> headers;
	
	while (i < argc) {
		std::string arg = argv[i];
		if (arg == "-f" || arg == "--file-format") {
			if (i+1 == argc) {
				std::cerr << "Option -f/--file-format needs a value!\n";
				return -1;
			}
			std::string val = argv[i+1];
			if (val == "plain" || val == "PLAIN") {
				file_format = FILE_FORMAT_PLAIN;
			} else if (val == "gzip" || val == "GZIP") {
				file_format = FILE_FORMAT_GZIP;
			} else if (val == "bin" || val == "BIN") {
				file_format = FILE_FORMAT_BIN;
			} else {
				std::cerr << "File format \"" << val
				          << " not recognized!\n";
				return -1;
			}
			i += 2;
		} else if (arg == "-i" || arg == "--info") {
			dump_info = true;
			++i;
		} else if (arg == "-n" || arg == "--nframes") {
			if (i+1 == argc) {
				std::cerr << "Option -n/--nframes requires a value!\n";
				return -1;
			}
			++i;
			frames = read_vector<int>(argc, argv, i);
			std::cerr << "After reading frames, i = " << i;
			if (i == argc) {
				std::cerr << " which is argc.\n";
			} else {
				std::cerr << " and arg there is " << argv[i] << "\n";
			}
		} else if (arg == "-c" || arg == "--columns") {
			if (i+1 == argc) {
				std::cerr << "Option -c/--columns requires a value!\n";
				return -1;
			}
			++i;
			headers = read_vector<std::string>(argc, argv, i);
			std::cerr << "After reading column headers, i = " << i;
			if (i == argc) {
				std::cerr << " which is argc.\n";
			} else {
				std::cerr << " and arg there is " << argv[i] << "\n";
			}
		} else {
			std::cerr << "Arg \"" << arg << "\" not recognized!\n";
			return -1;
		}
	}

	std::cerr << "Column headers: { ";
	for (const std::string c : headers) {
		std::cerr << c << " ";
	}
	std::cerr << "}\n";

	if (dump_info) {
		print_dump_info(dump_file, file_format, headers);
		return 0;
	}

	if (frames.empty()) {
		std::cerr << "No frames given!\n";
		return 0;
	}

	using namespace lammps_tools::readers;
	using namespace lammps_tools::writers;
	std::unique_ptr<dump_reader> dr(make_dump_reader(dump_file,
	                                                 file_format,
	                                                 DUMP_FORMAT_LAMMPS));
	dr->set_column_headers(headers);

	if (file_format == FILE_FORMAT_BIN) {
		int cur_frame = 0;
		int j = 0;
		block_data b;

		dump_reader_lammps_bin *d =
			dynamic_cast<dump_reader_lammps_bin*>(dr.get());

		for (int next_frame : frames) {
			d->skip_to_block(next_frame, cur_frame);
			if (d->next_block(b)) {
				std::cerr << "Error reading block "
				          << next_frame << "!\n";
				return -1;
			}
			std::cerr << "Read block at t=" << b.tstep << "\n";
			block_to_lammps_dump(std::cout, b, file_format);
			cur_frame = next_frame;
		}

		
	} else {
		int n_frames = 0;
		block_data b;
		while (dr->next_block(b) == 0) {
			if (util::contains(frames, n_frames)) {
				block_to_lammps_dump(std::cout, b, file_format);
			}
			++n_frames;
		}
	}
		
	
	return 0;
}
