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
  cpp_lib/bond_order.cpp
  cpp_lib/center_of_mass.cpp
  cpp_lib/correlation.cpp
  cpp_lib/cluster_finder.cpp
  cpp_lib/coarsening.cpp
  cpp_lib/data_field.cpp
  cpp_lib/data_reader_lammps.cpp
  cpp_lib/density_distribution.cpp
  cpp_lib/domain.cpp
  cpp_lib/dump_reader.cpp
  cpp_lib/dump_reader_hoomd_gsd.cpp
  cpp_lib/dump_reader_lammps.cpp
  cpp_lib/dump_reader_lammps_bin.cpp
  cpp_lib/dump_reader_lammps_gzip.cpp
  cpp_lib/dump_reader_lammps_plain.cpp
  cpp_lib/dump_reader_xyz.cpp
  cpp_lib/icosahedra.cpp
  cpp_lib/histogram.cpp
  cpp_lib/markov_state_capsid.cpp
  cpp_lib/msd.cpp
  cpp_lib/neighborize_bin.cpp
  cpp_lib/neighborize.cpp
  cpp_lib/neighborize_nsq.cpp
  cpp_lib/random_generator.cpp
  cpp_lib/rdf.cpp
  cpp_lib/scatter.cpp
  cpp_lib/skeletonize.cpp
  cpp_lib/topology.cpp
  cpp_lib/transformations.cpp
  cpp_lib/triangulate.cpp
  cpp_lib/util.cpp
  cpp_lib/writers_lammps.cpp)


set(CMAKE_SHARED_LINKER_FLAGS "-shared")

# Custom options:
option(USE_CGAL              "Use CGAL for Voronoi tesselation."            OFF)
option(USE_ARMADILLO         "Use Armadillo for normal mode analysis."      OFF)
option(USE_GSD               "Use GSD for reading HOOMD-Blue GSD files."    OFF)
option(USE_BOOST_GZIP        "Use boost for reading in GZIP files."         OFF)
option(USE_EXCEPTIONS        "Use C++ exceptions for error handling."       ON)
option(LEGACY_COMPILER       "Disable some features for ancient compilers." OFF)
option(USE_ASSERTIONS        "Compile library with assertions enabled."     ON)
option(THREADED_READ_BLOCKS  "Read blocks into a threaded queue."           OFF)

option(INCLUDE_C_INTERFACE "Compile the C interface into the library." ON)

if(ENABLE_MY_ASSERT)
  add_definitions(-DLT_DEBUG=1)
else()
  add_definitions(-DLT_DEBUG=0)
endif()

if(INCLUDE_C_INTERFACE)
  target_sources(lammpstools PRIVATE
    c_interface/lt_block_data.cpp
    c_interface/lt_block_writers.cpp
    c_interface/lt_data_field.cpp
    c_interface/lt_dump_reader.cpp
    c_interface/lt_transformations.cpp)
endif()


# Compiler settings:
set(CXX_STANDARD_REQUIRED "11")
set(CMAKE_CXX_NON_MS_FLAGS "-Werror=implicit -Werror=return-type -Werror=uninitialized -Wall")
set(CMAKE_CXX_FLAGS "-std=c++11 -g -O3 -pedantic -fPIC -ffast-math")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wno-missing-braces ${CMAKE_CXX_NON_MS_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -frounding-math ${CMAKE_CXX_NON_MS_FLAGS}")
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
    include_directories(${ARMADILLO_INCLUDE_DIR})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DHAVE_ARMADILLO")
  else()
    message(WARNING "Cannot find Armadillo! Not building with it!")
  endif()
endif(USE_ARMADILLO)


if(USE_GSD)
  set(GSD_DIR "../dependencies/gsd/")
  set(GSD_C_SOURCE "${GSD_DIR}/gsd.c" "../cpp_lib/writers_hoomd.cpp")
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



if(THREADED_READ_BLOCKS)
  set(READERWRITEQUEUE_DIR "../dependencies/readerwritequeue/")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTHREADED_READ_BLOCKS -I${READERWRITEQUEUE_DIR}")
endif(THREADED_READ_BLOCKS)


if(LEGACY_COMPILER)
  set(USE_EXCEPTIONS off)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLEGACY_COMPILER")
endif(LEGACY_COMPILER)

if(USE_EXCEPTIONS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_EXCEPTIONS")
endif(USE_EXCEPTIONS)


message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message("target_sources are ${target_sources}")
