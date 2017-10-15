#ifndef TYPES_HPP
#define TYPES_HPP

#include <algorithm>

#include <cstdint>

/*!
  \file  types.h
  \brief This file contains type definitions.

  \ingroup cpp_lib
*/

namespace lammps_tools {

typedef int64_t       bigint;    ///< signed long of guaranteed size (64 bits)

/// Use a union to easily read specific data from binary:
template <typename data_type>
union data_char_templ {
	data_type d;
	char c[sizeof(data_type)];
};


typedef data_char_templ<int32_t> utf32_char;
typedef data_char_templ<int16_t> utf16_char;
typedef data_char_templ<int8_t> utf8_char;
typedef data_char_templ<char> ascii_char;



}

#endif // TYPES_HPP
