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

cmake_minimum_required(VERSION 3.14) 

project(pteros)

# Optional dependencies
OPTION(WITH_PYTHON "Build python bindings and compile analysis plugins as Python modules" ON)
OPTION(WITH_OPENMP "Use OpenMP parallelization" ON)
OPTION(WITH_OPENBABEL "Use OpenBabel. Required to read pdbqt files and for substructure search." ON)
OPTION(WITH_GROMACS "Use Gromacs. Required to read tpr files." ON)
OPTION(WITH_TNG "Use TNG_IO. Required to read tng files." ON)

# Options to search for pre-installed dependencies
OPTION(TRY_SYSTEM_DEPENDENCIES "Try using system dependensies." ON)
OPTION(TRY_SYSTEM_EIGEN "Try using system Eigen." ON)
OPTION(TRY_SYSTEM_FMT "Try using system fmt." ON)
OPTION(TRY_SYSTEM_SCN "Try using system scn." ON)
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

#---------------------------------------------------

# Default to release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()

#-----------------------------
# Global compilation options:
#-----------------------------

# -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# C++17 and C99
set(CMAKE_CXX_STANDARD  17)
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

# Install headers
install(DIRECTORY include DESTINATION .)

# generate export targets file
install(
    EXPORT pterosTargets
    FILE pterosTargets.cmake
    NAMESPACE pteros::
    DESTINATION lib/cmake/pteros
)

# Generate cmake version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/pterosConfigVersion.cmake"
    VERSION 3.0
    COMPATIBILITY SameMajorVersion
)

# Generate cmake config file
SET(_plugin_dir ${PLUGINS_ABS_PATH})
configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/pterosConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/pterosConfig.cmake @ONLY
)

# Install CMake target machinery
install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/pterosConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/pterosConfigVersion.cmake
    DESTINATION
        lib/cmake/pteros
)

#---------------------------------------------
# Generate and install pterosrc file
#---------------------------------------------

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(_ld "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib64")
    SET(_py ${CMAKE_INSTALL_PREFIX}/python)
    SET(_path ${CMAKE_INSTALL_PREFIX}/bin)
    SET(_cmake ${CMAKE_INSTALL_PREFIX}/lib/cmake)
    configure_file(
        ${PROJECT_SOURCE_DIR}/cmake/pterosrc.in
        ${CMAKE_CURRENT_BINARY_DIR}/pterosrc @ONLY
    )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pterosrc DESTINATION lib)

    message("------------------------------")
    message("IMPORTANT!")
    message("In order to be able to link with` Pteros and to use Pteros Python modules")
    message("you have to source the pterosrc file in your bash shell:")
    message("source ${CMAKE_INSTALL_PREFIX}/lib/pterosrc")
    message("------------------------------")
ENDIF()
