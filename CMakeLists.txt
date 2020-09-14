#
# This file is a part of
#
# ============================================
#      Pteros molecular modeling library
# ============================================
#
# (C) 2009-2019, Semen Yesylevskyy
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

cmake_minimum_required(VERSION 3.1)

project(pteros)

OPTION(WITH_PYTHON "Build python bindings and compile analysis plugins as Python modules" ON)
OPTION(WITH_OPENMP "Use OpenMP parallelization" ON)
OPTION(WITH_OPENBABEL "Use OpenBabel. Required to read pdbqt files and for substructure search." ON)
OPTION(WITH_GROMACS "Use Gromacs. Required to read tpr files." ON)
OPTION(WITH_TNGIO "Use TNG_IO. Required to read tng files." ON)
OPTION(WITH_POWERSASA "Use POWERSASA code. This implies license restrictions described here: thirdparty/sasa/LICENSE" ON)

OPTION(MAKE_STANDALONE_PLUGINS "Compile analysis plugins as stand-alone executables" OFF)
OPTION(MAKE_PACKAGE "Generate package" OFF)
OPTION(MAKE_EXAMPLES "Compile examples and plugin templates" OFF)
OPTION(MAKE_TEST "Compile tests" OFF)

OPTION(DOWNLOAD_DEPS "Automatically download and compile dependencies if they are not found in the system" ON)

set(GROMACS_LIBRARIES "" CACHE STRING "Path to Gromacs libraries")
set(GROMACS_SOURCES "" CACHE STRING "Path to Gromacs source tree")

##############################
# Global compilation options:
##############################

# -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# C++14 and C99
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 99)

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

# Main Pteros include dir should be visible everywhere
include_directories(${PROJECT_SOURCE_DIR}/include)

# For downloading we need git
if(DOWNLOAD_DEPS)
    find_package(Git)
    if(NOT (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git"))
        message(INFO "Git is not found! Disabling automatic downloading of dependencies.")
        SET(DOWNLOAD_DEPS OFF)
    endif()
endif()

##############################
# Unconditional dependencies:
##############################

# Macro to download dependencies
MACRO(download module)
    if(DOWNLOAD_DEPS)
        message(STATUS "=======================")
        message(STATUS "Fetching git submodule for ${module}...")
        message(STATUS " ")
        execute_process(COMMAND git submodule update --init -- external/${module}
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
        if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/${module}/CMakeLists.txt")
            message(FATAL_ERROR "${module} was not downloaded!")
        endif()
        add_subdirectory(${PROJECT_SOURCE_DIR}/external/${module})
        message(STATUS "=======================")
    else()
        message(FATAL_ERROR "${module} was not found in the system. Try again with DOWNLOAD_DEPS=ON.")
    endif()
endmacro()


# Boost
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost 1.50 REQUIRED COMPONENTS system date_time filesystem)

# Eigen
find_package(Eigen3 3.3 NO_MODULE)
if(NOT Eigen3_FOUND)
    download(eigen)
    # Set Eigen location manually for openbabel.
    # By default it points to build dir which is empty until build stage
    set(Eigen3_location_for_babel ${PROJECT_SOURCE_DIR}/external/eigen)
endif()

# spdlog
find_package(spdlog CONFIG)
if(NOT spdlog_FOUND)
    set(SPDLOG_MASTER_PROJECT ON)
    set(SPDLOG_INSTALL ON)
    set(SPDLOG_BUILD_TESTS OFF)
    set(SPDLOG_BUILD_EXAMPLE OFF)
    download(spdlog)
endif()


##############################
# Conditional dependencies:
##############################

# OpenMP
if(WITH_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

# Python
if(WITH_PYTHON)    
    # Configure pybind11
    find_package(pybind11)
    if(NOT pybind11_FOUND)
        download(pybind11)
    endif()

    if(NOT MAKE_PACKAGE)
        # Set python install dir
        set(PY_INST_DIR python)
        set(PLUGINS_ABS_PATH ${CMAKE_INSTALL_PREFIX}/python/pteros_analysis_plugins)
    else()
        set(PY_INST_DIR ${PYTHON_SITE_PACKAGES})
        set(PLUGINS_ABS_PATH ${PY_INST_DIR}/pteros_analysis_plugins)
        if(DEFINED CMAKE_INSTALL_PREFIX)
            message(WARNING " You are building the package but set an install prefix.\n"
                            " Install prefix ${CMAKE_INSTALL_PREFIX} will be ignored.\n"
                            " Library files would be packaged to /usr\n"
                            " Python files would be packaged to ${PY_INST_DIR}\n"
                            " If this is not what you want rerun with -DMAKE_PACKAGE=OFF")
        endif()
    endif()
endif()

# OpenBabel
if(WITH_OPENBABEL)
    set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
    # Try to find OpenBabel 3
    find_package(OpenBabel3 3.0.0 REQUIRED)
    if(NOT OPENBABEL3_FOUND)        
        # Try to find OpenBabel 2
        message("OpenBabel v3 not found, searching for v2...")
        find_package(OpenBabel2 2.4.9 REQUIRED)
    endif()
endif()

# Gromacs
if(WITH_GROMACS)
    # See if pathes are provided
    if(NOT (GROMACS_SOURCES AND GROMACS_LIBRARIES))
        message(FATAL_ERROR "In order to link to Gromacs you need to provide path "
            "to Gromacs source tree matching your installed version "
            "using -DGROMACS_SOURCES=<path> "
            "and the path to installed Gromacs libs "
            "using -DGROMACS_LIBRARIES=<path> ")
    else()
        # Report what is used
        message("Searching for Gromacs sources in ${GROMACS_SOURCES}")
        message("Searching for Gromacs libs in ${GROMACS_LIBRARIES}")
        link_directories(${GROMACS_LIBRARIES})
    endif()
endif()


#----------------------------
# Build bundled thirdparty libraries
#----------------------------

# tng_io
if(WITH_TNGIO)
    add_subdirectory(thirdparty/tng_io)
endif()

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

#---------------------------------------------
# Compiling pteros itself
#---------------------------------------------

add_subdirectory(src/core)
add_subdirectory(src/analysis)
add_subdirectory(src/tools)
add_subdirectory(src/python)
add_subdirectory(src/python/compiled_analysis_plugins)
add_subdirectory(src/extras)

IF(MAKE_TEST)
    add_subdirectory(src/test)
ENDIF()

IF(MAKE_EXAMPLES)
    add_subdirectory(src/examples)
    add_subdirectory(template_plugin)
ENDIF()

#---------------------------------------------
# Installing
#---------------------------------------------

# Install headers
install(DIRECTORY include DESTINATION .)

# generate export targets file
install(EXPORT pterosTargets
    FILE pterosTargets.cmake
    NAMESPACE pteros::
    DESTINATION lib/cmake/pteros
)

# Generate cmake version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/cmake/pterosConfigVersion.cmake"
    VERSION 2.9
    COMPATIBILITY SameMajorVersion
)

# Generate cmake config file
SET(_plugin_dir ${PLUGINS_ABS_PATH})
configure_file(${PROJECT_SOURCE_DIR}/cmake/pterosConfig.cmake.in
               ${CMAKE_CURRENT_BINARY_DIR}/cmake/pterosConfig.cmake @ONLY)

# Install CMake target machinery
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/pterosConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/pterosConfigVersion.cmake
    DESTINATION
        lib/cmake/pteros)

#---------------------------------------------
# Generate and install pterosrc file
#---------------------------------------------

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(_ld "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib64")
    SET(_py ${CMAKE_INSTALL_PREFIX}/python)
    SET(_path ${CMAKE_INSTALL_PREFIX}/bin)
    configure_file(${PROJECT_SOURCE_DIR}/cmake/pterosrc.in
                   ${CMAKE_CURRENT_BINARY_DIR}/pterosrc @ONLY)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pterosrc DESTINATION lib)

    message("------------------------------")
    message("IMPORTANT!")
    message("In order to be able to link with` Pteros and to use Pteros Python modules")
    message("you have to source the pterosrc file in your bash shell:")
    message("source ${CMAKE_INSTALL_PREFIX}/lib/pterosrc")
    message("------------------------------")
ENDIF()

#--------------
# CPack config
#--------------
IF(MAKE_PACKAGE)
    message("Configuring package...")
    INCLUDE(InstallRequiredSystemLibraries)
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Pteros molecular modeling library")
    SET(CPACK_PACKAGE_VERSION_MAJOR "2")
    SET(CPACK_PACKAGE_VERSION_MINOR "9")
    SET(CPACK_PACKAGE_VERSION_PATCH "1")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        SET(CPACK_GENERATOR DEB)
        SET(CPACK_PACKAGE_CONTACT "Semen Yesylevskyy yesint4@yahoo.com")
        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "g++ (>=5.0), libboost-all-dev (>=1.50), python3-dev (>=3.5), python3-numpy, libeigen3-dev, libspdlog-dev, pybind11-dev")
    ENDIF()
    INCLUDE(CPack)
ENDIF()
