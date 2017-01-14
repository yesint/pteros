##########################################
# This file is MODIFIED by Pteros!
##########################################

message("Configuring TNG library...")

cmake_policy(SET CMP0022 OLD)

#------------------------------------------

if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
elseif(WIN32)
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /W2")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

option(BUILD_SHARED_LIBS "Enable shared libraries" ON)

option(TNG_BUILD_FORTRAN "Build Fortran compatible library and examples for testing" OFF)

option(TNG_BUILD_EXAMPLES "Build examples showing usage of the TNG API" ON)
option(TNG_BUILD_TEST "Build TNG testing binary." ON)
option(TNG_BUILD_COMPRESSION_TESTS "Build tests of the TNG compression library" OFF)

option(TNG_BUILD_OWN_ZLIB "Build and use the internal zlib library" OFF)
if(NOT TNG_BUILD_OWN_ZLIB)
  find_package(ZLIB QUIET)
endif()

include(CheckIncludeFile)
check_include_file(inttypes.h   HAVE_INTTYPES_H)

include(BuildTNG.cmake)
if (ZLIB_FOUND AND NOT TNG_BUILD_OWN_ZLIB)
  add_tng_io_library(tng_io)
else()
  add_tng_io_library(tng_io OWN_ZLIB)
endif()

# Use GNUInstallDirs to set paths on multiarch systems
include(GNUInstallDirs)


##########################################
# Install to lib directory of Pteros tree
##########################################

install(TARGETS tng_io
    RUNTIME DESTINATION lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
