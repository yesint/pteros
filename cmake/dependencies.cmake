#
# This file is a part of
#
# ============================================
#      Pteros molecular modeling library
# ============================================
#
# (C) 2009-2021, Semen Yesylevskyy
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

include(ExternalProject)
include(FetchContent)
cmake_policy(SET CMP0077 NEW) # To silence warnings

if(NOT TRY_SYSTEM_DEPENDENCIES)    
    set(TRY_SYSTEM_EIGEN       OFF)
    set(TRY_SYSTEM_SPDLOG      OFF)
    set(TRY_SYSTEM_FMT         OFF)
    set(TRY_SYSTEM_PYBIND11    OFF)
    set(TRY_SYSTEM_GROMACS     OFF)
    set(TRY_SYSTEM_OPENBABEL   OFF)
    message(STATUS "Searching for system-wide dependencies disbled.")
endif()

# To avoid updates on configure step
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
# To make it less noisy
set(FETCHCONTENT_QUIET ON)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# OpenMP
if(WITH_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
endif()

#======================================================
# Dependencies which use normal FetchContent workflow
#======================================================

# List of modules to fetch
set(fetch_list "")

#--------------------
# Eigen
#--------------------
if(TRY_SYSTEM_EIGEN)
    find_package(Eigen3 3.3 NO_MODULE)
endif()

if(NOT Eigen3_FOUND)
    if(NOT DOWNLOAD_DEPENDENCIES)
        message(FATAL_ERROR "Eigen3 is not available!")
    endif()

    message(STATUS "Will download Eigen")
    FetchContent_Declare(
      Eigen
      GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
      GIT_TAG           master
      GIT_SHALLOW       TRUE
      GIT_PROGRESS      TRUE
    )
    set(EIGEN_BUILD_DOC OFF CACHE INTERNAL "")
    set(BUILD_TESTING OFF CACHE INTERNAL "")
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE INTERNAL "")
    list(APPEND fetch_list Eigen)
endif()


#--------------------
# fmt
#--------------------
if(TRY_SYSTEM_FMT)
    find_package(fmt)
endif()

if(NOT fmt_FOUND)
    if(NOT DOWNLOAD_DEPENDENCIES)
        message(FATAL_ERROR "fmt is not available!")
    endif()

    message(STATUS "Will download and compile fmt in place")
    FetchContent_Declare(
            fmt
            GIT_REPOSITORY  https://github.com/fmtlib/fmt.git
            GIT_TAG         8.0.1
            GIT_SHALLOW     TRUE
            GIT_PROGRESS    TRUE
    )
    #set(FMT_MASTER_PROJECT ON CACHE INTERNAL "")
    set(FMT_INSTALL ON CACHE INTERNAL "")
    set(FMT_DOC OFF CACHE INTERNAL "")
    set(FMT_TEST OFF CACHE INTERNAL "")
    list(APPEND fetch_list fmt)
endif()


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
    FetchContent_Declare(
            spdlog
            GIT_REPOSITORY  https://github.com/gabime/spdlog.git
            GIT_TAG         v1.x
            GIT_SHALLOW     TRUE
            GIT_PROGRESS    TRUE
    )
    #set(SPDLOG_MASTER_PROJECT ON CACHE INTERNAL "")
    set(SPDLOG_INSTALL ON CACHE INTERNAL "")
    set(SPDLOG_BUILD_TESTS OFF CACHE INTERNAL "")
    set(SPDLOG_BUILD_EXAMPLE OFF CACHE INTERNAL "")
    set(SPDLOG_FMT_EXTERNAL ON CACHE INTERNAL "")
    list(APPEND fetch_list spdlog)    
endif()


#--------------------
# Python
#--------------------
if(WITH_PYTHON)    
    # Force Python3
    set(PYBIND11_PYTHON_VERSION 3 CACHE INTERNAL "")
    # Configure pybind11
    if(TRY_SYSTEM_PYBIND11)
        find_package(pybind11 QUIET)
    endif()

    if(NOT pybind11_FOUND)
        if(NOT DOWNLOAD_DEPENDENCIES)
            message(FATAL_ERROR "pybind11 is not available!")
        endif()

        message(STATUS "Will download and compile pybind11 in place")
        FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG        v2.2.3
            GIT_SHALLOW    TRUE
            GIT_PROGRESS   TRUE
        )        
        list(APPEND fetch_list pybind11)
    endif()

    # Set python install dir
    set(PY_INST_DIR python)
    set(PLUGINS_ABS_PATH ${CMAKE_INSTALL_PREFIX}/python/pteros_analysis_plugins)
endif()

# Fetch everything we need
if(fetch_list)
    message(STATUS "Will fetch the following: ${fetch_list}")
    FetchContent_MakeAvailable(${fetch_list})
endif()

#-----------------------------------------------------------
# Dependencies which use FetchContent for downloading only
# and ExternalProject for compiling with needed options
#-----------------------------------------------------------

# OpenBabel
include(${PROJECT_SOURCE_DIR}/cmake/openbabel.cmake)

# Gromacs
include(${PROJECT_SOURCE_DIR}/cmake/gromacs.cmake)

# TNG_IO
include(${PROJECT_SOURCE_DIR}/cmake/tngio.cmake)
