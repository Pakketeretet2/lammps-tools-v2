#include "lt_dump_reader.h"
#include "lt_block_data.h"

#include <vector>

int main( int argc, char **argv )
{
	lt_block_data_handle b;
	const char *fname = "A_32.0_R_40_ribbons.dump.bin";
	lt_dump_reader_handle handle = lt_new_dump_reader( fname, 2, 0 );
	std::vector<std::string> headers = { "id", "type", "x", "y", "z", "edt" };
	for( int i = 0; i < 6; ++i ){
		lt_set_col_header( handle, i, headers[i].c_str() );
	}

	for( int i = 0; i <= 20; ++i ){
		int status = lt_get_next_block( handle, &b );
		std::cerr << "i = " << i << ", status = " << status << ".\n";
	}

	return 0;
}
