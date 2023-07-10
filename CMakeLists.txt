#
# This file is a part of
#
# ============================================
#      Pteros molecular modeling library
# ============================================
#
# (C) 2009-2023, Semen Yesylevskyy
#
# All works, which use Pteros, should cite the following papers:
#
# 1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
#     molecular analysis library for C++ and python",
#     Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
#     doi: 10.1002/jcc.23943.
#
# 2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
#     library for molecular analysis",
#     Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
#     doi: 10.1002/jcc.22989.
#
# This is free software distributed under Artistic License:
# http://www.opensource.org/licenses/artistic-license-2.0.php
#

cmake_minimum_required(VERSION 3.18)

project(pteros VERSION 3.3 LANGUAGES CXX C)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR "In-source builds not allowed. Please make a build directory and run CMake from there.")
endif()

OPTION(MINIMAL_BUILD "Minimal build without heavy dependencies" OFF)

# Optional dependencies
OPTION(WITH_PYTHON "Build python bindings and compile analysis plugins as Python modules" ON)
OPTION(WITH_OPENMP "Use OpenMP parallelization" ON)

OPTION(WITH_OPENBABEL "Use OpenBabel. Required to read pdbqt files and for substructure search." ON)
OPTION(WITH_GROMACS "Use Gromacs. Required to read tpr files." ON)
OPTION(WITH_TNG "Use TNG_IO. Required to read tng files." ON)

# Options to search for pre-installed dependencies
OPTION(TRY_SYSTEM_DEPENDENCIES "Try using system dependensies." OFF)
OPTION(TRY_SYSTEM_EIGEN "Try using system Eigen." ON)
OPTION(TRY_SYSTEM_FMT "Try using system fmt." ON)
OPTION(TRY_SYSTEM_SPDLOG "Try using system spdlog." ON)
OPTION(TRY_SYSTEM_PYBIND11 "Try using system pybind11." ON)
OPTION(TRY_SYSTEM_GROMACS "Try using system Gromacs." ON)
OPTION(TRY_SYSTEM_OPENBABEL "Try using system OpenBabel" ON)

OPTION(DOWNLOAD_DEPENDENCIES "Download and compile dependencies if not found in the system" ON)

OPTION(MAKE_TOOLS "Build pteros tools" ON)
OPTION(MAKE_STANDALONE_PLUGINS "Compile analysis plugins as stand-alone executables" OFF)
OPTION(MAKE_PACKAGE "Generate package" OFF)
OPTION(MAKE_EXAMPLES "Compile examples and plugin templates" OFF)
OPTION(MAKE_TEST "Compile tests" OFF)

if(MINIMAL_BUILD)
    message(STATUS "This is a minimal build without Gromacs, OpenBabel and TNG")
    set(WITH_GROMACS OFF CACHE BOOL "" FORCE)
    set(WITH_OPENBABEL OFF CACHE BOOL "" FORCE)
    set(WITH_TNG OFF CACHE BOOL "" FORCE)
endif()

#---------------------------------------------------

# Default to release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()

# Manage installation location
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # No custom install prefix is given
    # Check if we are inside the venv, if so install to this env
    if(DEFINED ENV{VIRTUAL_ENV})
        set(CMAKE_INSTALL_PREFIX $ENV{VIRTUAL_ENV} CACHE PATH "Choose install location" FORCE)
    else()
        if(UNIX AND NOT APPLE)
            set(CMAKE_INSTALL_PREFIX "~/.local" CACHE PATH "Choose install location" FORCE)
        endif()
    endif()
endif()

#-----------------------------
# Global compilation options:
#-----------------------------

# Use ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# C++17 and C99
set(CMAKE_CXX_STANDARD  20)
set(CMAKE_C_STANDARD    99)

#Sets optmization in Release build
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

# Set specific options for GCC if used
if(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
endif()

# Strip libraries and executables to decrease size in Release builds
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -s")
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s")
ENDIF()

#-------------------------------------------------------------------------------
# Interface library to provide basic includes for other targets at build time
#-------------------------------------------------------------------------------

add_library(pteros_interface INTERFACE)
target_include_directories(pteros_interface INTERFACE ${PROJECT_SOURCE_DIR}/include/)

#----------------------------
# Manage dependencies
#----------------------------

# Set CPM cache location if not provide in the command line
if(NOT CPM_SOURCE_CACHE)
    set(CPM_SOURCE_CACHE "~/.local/CPM_CACHE" CACHE STRING "Location of CPM cache" FORCE)
endif()

# include CPM
include(${PROJECT_SOURCE_DIR}/cmake/CPM.cmake)

# Get modified PackageProject.cmake
CPMAddPackage("gh:yesint/PackageProject.cmake#master")

# Set versions of dependencies
set(EIGEN_VERSION 3.4)
set(FMT_VERSION 8.1.1)
set(SPDLOG_VERSION 1.12.0)

# Manage dependencies
include(${PROJECT_SOURCE_DIR}/cmake/dependencies.cmake)

#------------------------------------
# Build bundled thirdparty libraries
#------------------------------------

# DSSP library
add_subdirectory(thirdparty/dssp)
# VMD molfile plugins
add_subdirectory(thirdparty/molfile_plugins)
# xdrfile library
add_subdirectory(thirdparty/xdrfile)
# voro++ library
add_subdirectory(thirdparty/voro++)
# SASA library from MDTraj
add_subdirectory(thirdparty/sasa)

#----------------------------
# Compiling pteros itself
#----------------------------
add_subdirectory(src/core)
add_subdirectory(src/analysis)
add_subdirectory(src/python)
add_subdirectory(src/python/compiled_analysis_plugins)
add_subdirectory(src/extras)

IF(MAKE_TOOLS)
    add_subdirectory(src/tools)
ENDIF()

IF(MAKE_TEST)
    add_subdirectory(src/test)
ENDIF()

IF(MAKE_EXAMPLES)
    add_subdirectory(src/examples)
    add_subdirectory(template_plugin)
ENDIF()

#---------------------
# Installing
#---------------------
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
if(WITH_PYTHON)
    message(STATUS "Python install dir: ${PY_INST_DIR}")
endif()

packageProject(
  # the name of the target to export
  NAME ${PROJECT_NAME}
  # the version of the target to export
  VERSION ${PROJECT_VERSION}
  # a temporary directory to create the config files
  BINARY_DIR ${PROJECT_BINARY_DIR}
  # location of the target's public headers
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}
  # should match the target's INSTALL_INTERFACE include directory
  INCLUDE_DESTINATION include/${PROJECT_NAME}
  # semicolon separated list of the project's dependencies
  DEPENDENCIES "Eigen3 ${EIGEN_VERSION};fmt ${FMT_VERSION};spdlog ${SPDLOG_VERSION}"
  # (optional) create a header containing the version info
  # Note: that the path to headers should be lowercase
  VERSION_HEADER "${PROJECT_NAME}/version.h"
  # (optional) install your library with a namespace (Note: do NOT add extra '::')
  NAMESPACE pteros
  # (optional) define the project's version compatibility, defaults to `AnyNewerVersion`
  # supported values: `AnyNewerVersion|SameMajorVersion|SameMinorVersion|ExactVersion`
  COMPATIBILITY AnyNewerVersion
  # (optional) option to disable the versioning of install destinations
  DISABLE_VERSION_SUFFIX YES
  # Automatically install all targets with TO_INSTALL property set to YES
  AUTO_FIND_TARGETS_TO_INSTALL YES
)

#---------------------------------------------
# Generate and install pterosrc file
#---------------------------------------------

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(_ld "${CMAKE_INSTALL_PREFIX}/lib/pteros:${CMAKE_INSTALL_PREFIX}/lib64")
    SET(_py ${PY_INST_DIR})
    SET(_path ${CMAKE_INSTALL_PREFIX}/bin/pteros)
    SET(_cmake ${CMAKE_INSTALL_PREFIX}/lib/cmake)
    configure_file(
        ${PROJECT_SOURCE_DIR}/cmake/pterosrc.in
        ${CMAKE_CURRENT_BINARY_DIR}/pterosrc @ONLY
    )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pterosrc DESTINATION lib/pteros)

    message("------------------------------")
    message("IMPORTANT!")
    message("In order to be able to link with` Pteros and to use Pteros Python modules")
    message("you have to source the pterosrc file in your bash shell:")
    message("source ${CMAKE_INSTALL_PREFIX}/lib/pteros/pterosrc")
    message("------------------------------")
ENDIF()
