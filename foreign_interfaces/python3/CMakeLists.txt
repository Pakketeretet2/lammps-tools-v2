cmake_minimum_required (VERSION 2.8)
project (lammps-tools-python-interface)
# The version number.
set(Tutorial_VERSION_MAJOR 2)
set(Tutorial_VERSION_MINOR 0)

# Compiler settings:
set(CXX_STANDARD_REQUIRED "11")
set(CMAKE_CXX_NON_MS_FLAGS "-Werror=implicit -Werror=return-type -Werror=uninitialized -Wall")
set(CMAKE_CXX_FLAGS "-std=c++11 -g -O3 -pedantic -fPIC")

include_directories( ../../c_interface/ ../../dependencies/ )


if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wno-missing-braces ${CMAKE_CXX_NON_MS_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -frounding-math ${CMAKE_CXX_NON_MS_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

endif()

if(UNIX)
  set(EXTENSION "_.so")
else()
  set(EXTENSION "_")
endif()

# Check if you find Python at all:
find_package( PythonLibs )
if( PythonLibs_FOUND )
  message("Found Python library: ${PYTHON_LIBRARIES}")
  message("Found Python header in ${PYTHON_INCLUDE_DIRS}")
  include_directories( AFTER ${PYTHON_INCLUDE_DIRS} )
else()
  message(SEND_WARNING " Could not automatically find Python! Set it yourself in advanced options!")
endif()

include_directories( AFTER ${PYTHON_INCLUDE_DIRS} )

set(LAMMPSTOOLS_DIR "../../" CACHE FILEPATH "Path to dir containing liblammpstools.so")
find_library(LAMMPSTOOLS_LIBRARY lammpstools HINTS ${LAMMPSTOOLS_DIR})
set(TARGET_LIBS ${LAMMPSTOOLS_LIBRARY})

set(PYBIND11_DIR "pybind11")
add_subdirectory( ${PYBIND11_DIR} )

execute_process( COMMAND "python-config --cflags --ldflags" OUTPUT_VARIABLE PYTHON_FLAGS )

foreach(SRC_FILE
    lammpstools/block_data_python.cpp
    lammpstools/block_writers_python.cpp
    lammpstools/bond_order_python.cpp
    lammpstools/conversions_python.cpp
    lammpstools/dump_reader_python.cpp
    lammpstools/data_field_python.cpp
    lammpstools/neighborize_python.cpp
    lammpstools/scatter_python.cpp
    lammpstools/transformations_python.cpp)

  string(REPLACE "lammpstools/" "" TMP ${SRC_FILE})
  string(REPLACE "_python.cpp" "${EXTENSION}" OBJ_NAME ${TMP})

  pybind11_add_module( ${OBJ_NAME} ${SRC_FILE} ${TARGET_LIBS} )
  target_link_libraries( ${OBJ_NAME} PRIVATE ${TARGET_LIBS} )

endforeach(SRC_FILE)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PYTHON_FLAGS}")
