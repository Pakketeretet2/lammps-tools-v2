#ifndef ENUMS_HPP
#define ENUMS_HPP

/**
   \file enums.hpp

   This file contains some globally-visible enums.
*/

namespace lammps_tools {

/// Enumerates atom styles
enum LT_ATOM_STYLES {
	ATOM_STYLE_ATOMIC = 0, ///< Atomic data only
	ATOM_STYLE_MOLECULAR   ///< Data on molecules too
};

/**
   \brief Enumerates the file formats.
*/
enum LT_DUMP_READER_FILE_FORMATS {
	FILE_FORMAT_UNSET = -1, ///< Not set/invalid format
	FILE_FORMAT_PLAIN = 0,  ///< Plain text
	FILE_FORMAT_GZIP,       ///< Gzipped plain text
	FILE_FORMAT_BIN         ///< Binary format

};

/**
   \brief Enumerates dump formats.
*/
enum LT_DUMP_READER_DUMP_FORMATS {
	DUMP_FORMAT_UNSET = -1,    ///< Not set/invalid format
	DUMP_FORMAT_LAMMPS = 0,    ///< LAMMPS dump file
	DUMP_FORMAT_LAMMPS_LOCAL,  ///< LAMMPS dump file
	DUMP_FORMAT_HOOMD,         ///< HOOMD-blue GSD dump file
	DUMP_FORMAT_NAMD,          ///< NAMD style DCD format
	DUMP_FORMAT_XYZ            ///< Standard XYZ dump file.
};

}

#endif // ENUMS_HPP
