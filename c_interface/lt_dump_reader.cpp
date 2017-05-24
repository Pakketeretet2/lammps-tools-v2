#include "lt_dump_reader.h"

#include "../cpp_lib/enums.hpp"
#include "../cpp_lib/dump_reader_lammps.hpp"
#include "../cpp_lib/block_data_access.hpp"


// ****** Helper functions:  ********
lammps_tools::readers::dump_reader_lammps
*attempt_lammps_dump_reader_cast( lt_dump_reader_handle drh )
{
	using lammps_tools::readers::dump_reader_lammps;
	dump_reader_lammps *drl = static_cast<dump_reader_lammps*>( drh.dr );
	if( !drl ){
		std::cerr << "Error casting to LAMMPS dump reader!\n";
		return nullptr;
	}
	return drl;
}


// ******  Interface implementation:  ********



extern "C" {

lt_dump_reader_handle lt_new_dump_reader( const char *fname,
                                          int fformat, int dformat )
{
	lt_dump_reader_handle drh;
	drh.dr = lammps_tools::readers::make_dump_reader( fname, fformat, dformat );
	drh.fformat = fformat;
	drh.dformat = dformat;
	std::cerr << "Created new dump_reader at " << drh.dr
	          << " for fformat = " << fformat << " and dformat = "
	          << dformat << ".\n";
	return drh;
}



lt_dump_reader_handle lt_new_dump_reader_local( const char *fname,
                                                int fformat, int dformat )
{
	lt_dump_reader_handle drh;
	int dstyle = lammps_tools::readers::dump_reader_lammps::LOCAL;

	if( fformat != lammps_tools::FILE_FORMAT_PLAIN ){
		drh.dr = nullptr;
		return drh;
	}

	drh.dr = lammps_tools::readers::make_dump_reader_lammps( fname, fformat,
	                                                         dstyle );
	drh.fformat = fformat;
	drh.dformat = dformat;
	std::cerr << "Created new dump_reader_lammps at " << drh.dr
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
	if( dr == nullptr ){
		return POINTER_NULL;
	}

	if( dr->eof() )  return AT_EOF;
	if( dr->good() ) return IS_GOOD;
	else             return IS_BAD;

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
	if( drh.dformat != lammps_tools::DUMP_FORMAT_LAMMPS ){
		std::cerr << "Ignoring column headers for non-LAMMPS dump...\n";
		return;
	}
	if( dump_reader_lammps *drl = attempt_lammps_dump_reader_cast( drh ) ){
		drl->set_column_header( n, header );
	}
}

bool lt_set_column_header_as_special( lt_dump_reader_handle drh,
                                      const std::string &header,
                                      int special_field_type )
{
	using lammps_tools::readers::dump_reader_lammps;
	if( drh.dformat != lammps_tools::DUMP_FORMAT_LAMMPS ){
		std::cerr << "Ignoring column headers for non-LAMMPS dump...\n";
		return false;
	}

	if( dump_reader_lammps *drl = attempt_lammps_dump_reader_cast( drh ) ){
		return drl->set_column_header_as_special( header,
		                                          special_field_type );
	}else{
		std::cerr << "Error casting to LAMMPS dump reader!\n";
		return false;
	}
}

bool lt_set_column_type( lt_dump_reader_handle drh,
                         const std::string &header, int type )
{
	using lammps_tools::readers::dump_reader_lammps;
	if( drh.dformat != lammps_tools::DUMP_FORMAT_LAMMPS ){
		std::cerr << "Ignoring column headers for non-LAMMPS dump...\n";
		return false;
	}

	if( dump_reader_lammps *drl = attempt_lammps_dump_reader_cast( drh ) ){
		drl->set_column_type( header, type );
		return true;
	}else{
		std::cerr << "Error casting to LAMMPS dump reader!\n";
		return false;
	}
}

int  lt_get_column_type( lt_dump_reader_handle drh,
                         const std::string &header )
{
	using lammps_tools::readers::dump_reader_lammps;
	if( drh.dformat != lammps_tools::DUMP_FORMAT_LAMMPS ){
		std::cerr << "Ignoring column headers for non-LAMMPS dump...\n";
		return -1;
	}

	if( dump_reader_lammps *drl = attempt_lammps_dump_reader_cast( drh ) ){
		return drl->get_column_type( header );
	}else{
		std::cerr << "Error casting to LAMMPS dump reader!\n";
		return -1;
	}
}

void lt_set_default_column_type( lt_dump_reader_handle drh, int type )
{
	using lammps_tools::readers::dump_reader_lammps;
	if( dump_reader_lammps *drl = attempt_lammps_dump_reader_cast( drh ) ){
		drl->set_default_column_type( type );
	}else{
		std::cerr << "Error casting to LAMMPS dump reader!\n";
	}

}

} // extern "C"
