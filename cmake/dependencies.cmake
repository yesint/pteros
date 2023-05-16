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
#---------------------------------------------------

cmake_policy(SET CMP0077 NEW) # To silence warnings

if(NOT TRY_SYSTEM_DEPENDENCIES)    
    set(TRY_SYSTEM_EIGEN       OFF)
    set(TRY_SYSTEM_SPDLOG      OFF)
    set(TRY_SYSTEM_FMT         OFF)
    #set(TRY_SYSTEM_SCN         OFF)
    set(TRY_SYSTEM_PYBIND11    OFF)
    set(TRY_SYSTEM_GROMACS     OFF)
    set(TRY_SYSTEM_OPENBABEL   OFF)
    message(STATUS "Searching for system-wide dependencies disbled.")
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# OpenMP
if(WITH_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        # Link OpenMP globally to all pteros targets
        link_libraries(OpenMP::OpenMP_CXX)
    else()
        message(WARNING "OpenMP is not available.")
    endif()
endif()

#======================================================
# Dependencies which can use CPM for installing
#======================================================

include(${PROJECT_SOURCE_DIR}/cmake/CPM.cmake)

#--------------------
# Eigen
#--------------------
if(TRY_SYSTEM_EIGEN)
    find_package(Eigen3 3.4 NO_MODULE)
endif()

if(NOT Eigen3_FOUND)
    if(NOT DOWNLOAD_DEPENDENCIES)
        message(FATAL_ERROR "Eigen3 is not available!")
    endif()

    message(STATUS "Will download Eigen")
    CPMAddPackage(
        NAME              Eigen
        GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
        GIT_TAG           3.4 #master
        OPTIONS
            "EIGEN_BUILD_DOC OFF"
            "BUILD_TESTING OFF"
            "EIGEN_BUILD_TESTING OFF"
            "EIGEN_BUILD_PKGCONFIG OFF"
    )
endif()


#--------------------
# fmt
#--------------------
if(TRY_SYSTEM_FMT)
    find_package(fmt QUIET)
endif()

if(NOT fmt_FOUND)
    if(NOT DOWNLOAD_DEPENDENCIES)
        message(FATAL_ERROR "fmt is not available!")
    endif()

    message(STATUS "Will download and compile fmt in place")
    CPMAddPackage(
        NAME                fmt
        GITHUB_REPOSITORY   fmtlib/fmt
        GIT_TAG             8.1.1
        OPTIONS
            "FMT_INSTALL ON"
            "FMT_DOC OFF"
            "FMT_TEST OFF"
    )
endif()

#--------------------
# scn
#--------------------
#if(TRY_SYSTEM_SCN)
#    find_package(scn QUIET)
#endif()
#
#if(NOT scn_FOUND)
#    if(NOT DOWNLOAD_DEPENDENCIES)
#        message(FATAL_ERROR "scn is not available!")
#    endif()
#
#    message(STATUS "Will download and compile scn in place")
#    CPMAddPackage(
#        NAME                scn
#        GITHUB_REPOSITORY   eliaskosunen/scnlib
#        GIT_TAG             v1.1.2
#        OPTIONS
#            "SCN_INSTALL ON"
#            "SCN_DOCS OFF"
#            "SCN_TESTS OFF"
#            "SCN_EXAMPLES OFF"
#            "SCN_BENCHMARKS OFF"
#    )
#endif()

#--------------------
# spdlog
#--------------------
if(TRY_SYSTEM_SPDLOG)
    find_package(spdlog CONFIG)
endif()

if(NOT spdlog_FOUND)
    if(NOT DOWNLOAD_DEPENDENCIES)
        message(FATAL_ERROR "spdlog is not available!")
    endif()

    message(STATUS "Will download and compile spdlog in place")
    CPMAddPackage(
        NAME                spdlog
        GITHUB_REPOSITORY   gabime/spdlog
        GIT_TAG             v1.x
        OPTIONS
            "SPDLOG_INSTALL ON"
            "SPDLOG_BUILD_TESTS OFF"
            "SPDLOG_BUILD_EXAMPLE OFF"
            "SPDLOG_FMT_EXTERNAL ON"
    )
endif()


#--------------------
# Pybind11
#--------------------
if(WITH_PYTHON)    
    # Force Python3
    set(PYBIND11_PYTHON_VERSION 3 CACHE INTERNAL "")
    # Configure pybind11

    # By default we will download it.
    # We have to call CPM each time when configuring without system pybind11
    # Even if pybind11_FOUND was already set previously!
    # Without this pybind11 command are not recognized
    # Only if found in system, this will be disabled
    set(DOWNLOAD_PYBIND ON)

    if(TRY_SYSTEM_PYBIND11)
        find_package(pybind11 QUIET)
        if(NOT pybind11_FOUND)
            if(NOT DOWNLOAD_DEPENDENCIES)
                message(FATAL_ERROR "pybind11 is not available!")
            endif()
        else()
            # We found system pybind, disable download
            set(DOWNLOAD_PYBIND OFF)
        endif()
    endif()

    # Always call CPM unless explicitly disabled
    if(DOWNLOAD_PYBIND)
        CPMAddPackage(
            NAME                pybind11
            GITHUB_REPOSITORY   pybind/pybind11
            GIT_TAG             master
        )
    endif()

    # Set python install dir
    set(PY_INST_DIR python)
    set(PLUGINS_ABS_PATH ${CMAKE_INSTALL_PREFIX}/python/pteros_analysis_plugins)
endif()


#-----------------------------------------------------------------
# Dependencies which use ExternalProject for compiling
#-----------------------------------------------------------------
include(ExternalProject)

# OpenBabel
include(${PROJECT_SOURCE_DIR}/cmake/openbabel.cmake)
# Gromacs
include(${PROJECT_SOURCE_DIR}/cmake/gromacs.cmake)
# TNG_IO
include(${PROJECT_SOURCE_DIR}/cmake/tngio.cmake)
