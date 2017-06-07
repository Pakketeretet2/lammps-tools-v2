cmake_minimum_required (VERSION 2.8)
project (lammps-tools)
# The version number.
set(lammpstools_VERSION_MAJOR 2)
set(lammpstools_VERSION_MINOR 0)

set(CXX_STANDARD_REQUIRED "11")

# Set up source files:
add_library(lammpstools SHARED
  cpp_lib/atom_type_info.cpp
  cpp_lib/block_data.cpp
  cpp_lib/cluster_finder.cpp
  cpp_lib/data_field.cpp
  cpp_lib/data_reader_lammps.cpp
  cpp_lib/domain.cpp
  cpp_lib/dump_reader.cpp
  cpp_lib/dump_reader_lammps.cpp
  cpp_lib/dump_reader_lammps_bin.cpp
  cpp_lib/dump_reader_lammps_gzip.cpp
  cpp_lib/dump_reader_lammps_plain.cpp
  cpp_lib/markov_state_model.cpp
  cpp_lib/neighborize_bin.cpp
  cpp_lib/neighborize.cpp
  cpp_lib/neighborize_nsq.cpp
  cpp_lib/rdf.cpp
  cpp_lib/skeletonize.cpp
  cpp_lib/topology.cpp
  cpp_lib/triangulate.cpp
  cpp_lib/util.cpp
  cpp_lib/writers_lammps.cpp
  c_interface/lt_block_data.cpp
  c_interface/lt_data_field.cpp
  c_interface/lt_dump_reader.cpp)

set(CMAKE_SHARED_LINKER_FLAGS "-shared")

# Custom options:
option(USE_CGAL        "Use CGAL for Voronoi tesselation."            OFF)
option(USE_ARMADILLO   "Use Armadillo for normal mode analysis."      OFF)
option(USE_GSD         "Use GSD for reading HOOMD-Blue GSD files."    OFF)
option(USE_BOOST_GZIP  "Use boost for reading in GZIP files."         OFF)
option(USE_EXCEPTIONS  "Use C++ exceptions for error handling."       ON)
option(LEGACY_COMPILER "Disable some features for ancient compilers." OFF)

# Compiler settings:
set(CMAKE_CXX_FLAGS "-std=c++11 -g -O3 -pedantic -fPIC -Wall -Werror=implicit")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=uninitialized -Weffc++")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wno-missing-braces")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -frounding-math")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

endif()


## Handle custom options:
if(USE_CGAL)
  find_package(CGAL)
  if(CGAL_FOUND)
    target_link_libraries(lammpstools -lCGAL -lgmp)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DHAVE_LIB_CGAL")
  else()
    message(WARNING "Cannot find CGAL! Not building with it!")
  endif()
endif(USE_CGAL)


if(USE_ARMADILLO)
  find_package(Armadillo)
  if(Armadillo_FOUND)
    target_link_libraries(lammpstools -larmadillo)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DHAVE_ARMADILLO")
  else()
    message(WARNING "Cannot find Armadillo! Not building with it!")
  endif()
endif(USE_ARMADILLO)


if(USE_GSD)
  set(GSD_DIR "dependencies/gsd/")
  set(GSD_C_SOURCE "${GSD_DIR}/gsd.c")
  target_sources(lammpstools PRIVATE ${GSD_C_SOURCE})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DHAVE_GSD -I${GSD_DIR}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
endif(USE_GSD)


if(USE_BOOST_GZIP)
  find_package(Boost)
  if(Boost_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DHAVE_BOOST_GZIP")
    target_link_libraries(lammpstools -lboost_iostreams)
  else()
    message(WARNING "Cannot find Boost! Not building with it!")
  endif()
endif(USE_BOOST_GZIP)

if(LEGACY_COMPILER)
  set(USE_EXCEPTIONS off)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLEGACY_COMPILER")
endif(LEGACY_COMPILER)

if(USE_EXCEPTIONS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_EXCEPTIONS")
endif(USE_EXCEPTIONS)
