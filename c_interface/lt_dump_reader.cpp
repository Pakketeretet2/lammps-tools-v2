#include "lt_dump_reader.h"

#include "../cpp_lib/dump_reader_lammps.hpp"

extern "C" {

lt_dump_reader_handle lt_new_dump_reader( const char *fname,
                                          int fformat, int dformat )
{
	lt_dump_reader_handle drh;
	drh.dr = lammps_tools::readers::make_dump_reader( fname,
	                                                  dformat, fformat );
	drh.fformat = fformat;
	drh.dformat = dformat;
	std::cerr << "Created new dump_reader at " << drh.dr
	          << " for fformat = " << fformat << " and dformat = "
	          << dformat << ".\n";
	return drh;
}

void lt_delete_dump_reader( lt_dump_reader_handle drh )
{
	delete drh.dr;
	std::cerr << "Deleted dump_reader at " << drh.dr << ".\n";
}


int lt_dump_reader_status( lt_dump_reader_handle drh )
{
	lammps_tools::readers::dump_reader *dr = drh.dr;
	if( dr == nullptr ) return POINTER_NULL;

	if     ( dr->good() ) return IS_GOOD;
	else if( dr->eof()  ) return AT_EOF;
	else                  return IS_BAD;
}

int lt_get_next_block( lt_dump_reader_handle drh, lt_block_data_handle *bdh )
{
	//lammps_tools::block_data tmp;
	//*bdh->bd = tmp;
	return drh.dr->next_block( *bdh->bd );
}



int lt_number_of_blocks( lt_dump_reader_handle drh )
{
	return number_of_blocks( *(drh.dr) );
}


void lt_set_col_header( lt_dump_reader_handle drh, int n, const char *header )
{
	using lammps_tools::readers::dump_reader_lammps;
	if( drh.dformat != lammps_tools::readers::LAMMPS ){
		std::cerr << "Ignoring column headers for non-LAMMPS dump...\n";
		return;
	}

	dump_reader_lammps *drl = static_cast<dump_reader_lammps*>( drh.dr );
	drl->set_column_header( n, header );
}

} // extern "C"
