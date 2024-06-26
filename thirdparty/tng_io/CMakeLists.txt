cmake_minimum_required(VERSION 3.1)

project(TNG_IO)

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

option(TNG_BUILD_EXAMPLES "Build examples showing usage of the TNG API" OFF)
option(TNG_BUILD_TEST "Build TNG testing binary." OFF)
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

# Create the tng_ioConfig.cmake and tng_ioConfigVersion.cmake files for the install tree
configure_file(              src/lib/tng_io-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/tng_io-config.cmake" @ONLY)
configure_file(              src/lib/tng_io-configVersion.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/tng_io-configVersion.cmake" @ONLY)

install(TARGETS tng_io
        EXPORT tng_io
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(EXPORT tng_io FILE tng_io.cmake
        NAMESPACE tng_io::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/tng_io")

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/tng_io-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/tng_io-configVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/tng_io")

install(FILES
        include/tng/tng_io.h include/tng/tng_io_fwd.h
        #${CMAKE_CURRENT_BINARY_DIR}/include/tng/version.h
        ${TNG_ROOT_BINARY_DIR}/include/tng/version.h
        DESTINATION include/tng)

#-- Add an Option to toggle the generation of the API documentation
option(TNG_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(TNG_BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
  #-- Configure the Template Doxyfile for our specific project
  configure_file(Doxyfile.in
                 ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
  #-- Add a custom target to run Doxygen when ever the project is built
  add_custom_target (Docs ALL
                                        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
                                        SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.

  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Documentation/
      DESTINATION share/tng/doc)
endif()

add_subdirectory(src)
