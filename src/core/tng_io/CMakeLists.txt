##########################################
# This file is MODIFIED by Pteros!
##########################################

message("Configuring TNG library...")

if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
elseif(WIN32 AND NOT CYGWIN)
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /W2")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

option(BUILD_SHARED_LIBS "Enable shared libraries" ON)

option(TNG_BUILD_FORTRAN "Build Fortran compatible library and examples for testing" OFF)

option(TNG_BUILD_EXAMPLES "Build examples showing usage of the TNG API" ON)
option(TNG_BUILD_TEST "Build TNG testing binary." ON)
option(TNG_BUILD_COMPRESSION_TESTS "Build tests of the TNG compression library" OFF)

find_package(ZLIB QUIET)
option(TNG_BUILD_WITH_ZLIB "Build TNG with zlib compression" ${ZLIB_FOUND})

include(CheckIncludeFile)
check_include_file(inttypes.h   HAVE_INTTYPES_H)

include(BuildTNG.cmake)
tng_get_source_list(TNG_SOURCES TNG_COMPILE_DEFS)

tng_set_source_properties(WITH_ZLIB ${ZLIB_FOUND})

if (ZLIB_FOUND)
  list(APPEND EXTRA_LIBRARIES ${ZLIB_LIBRARIES})
endif()

if (UNIX)
  list(APPEND EXTRA_LIBRARIES m)
endif()

add_library(tng_io ${TNG_SOURCES})

target_link_libraries(tng_io ${EXTRA_LIBRARIES})

set_target_properties(tng_io PROPERTIES VERSION ${TNG_IO_VERSION} SOVERSION ${TNG_MAJOR_VERSION})

##########################################
# Install to lib directory of Pteros tree
##########################################

install(TARGETS tng_io
    RUNTIME DESTINATION lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
