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
    set(TRY_SYSTEM_PYBIND11    OFF)
    set(TRY_SYSTEM_GROMACS     OFF)
    set(TRY_SYSTEM_OPENBABEL   OFF)
    message(STATUS "Searching for system-wide dependencies disbled.")
endif()

# To avoid updates on configure step
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
set(FETCHCONTENT_QUIET OFF)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

##############################
# Unconditional dependencies:
##############################

# List of modules to fetch
set(fetch_list "")

#--------------------
# Boost
#--------------------
#set(Boost_USE_STATIC_LIBS OFF)
#find_package(Boost 1.50 REQUIRED COMPONENTS system date_time filesystem)

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
    set(SPDLOG_MASTER_PROJECT ON CACHE INTERNAL "")
    set(SPDLOG_INSTALL ON CACHE INTERNAL "")
    set(SPDLOG_BUILD_TESTS OFF CACHE INTERNAL "")
    set(SPDLOG_BUILD_EXAMPLE OFF CACHE INTERNAL "")
    list(APPEND fetch_list spdlog)
    #FetchContent_MakeAvailable(spdlog)
endif()

##############################
# Conditional dependencies:
##############################

# OpenMP
if(WITH_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
endif()

#--------------------
# Python
#--------------------
if(WITH_PYTHON)    
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
        # Force Python3
        set(PYBIND11_PYTHON_VERSION 3 CACHE INTERNAL "")
        list(APPEND fetch_list pybind11)
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

# Fetch everything we need
if(fetch_list)
    message(STATUS "Will fetch the following: ${fetch_list}")
    FetchContent_MakeAvailable(${fetch_list})
endif()

#--------------------
# OpenBabel
#--------------------
if(WITH_OPENBABEL)
    if(TRY_SYSTEM_OPENBABEL)        
        # Try to find OpenBabel 3
        find_package(OpenBabel3 3.0.0)
        if(NOT OPENBABEL3_FOUND)
            # Try to find OpenBabel 2
            message(STATUS "OpenBabel v3 not found, searching for v2...")
            find_package(OpenBabel2 2.4.9)
            if(NOT OPENBABEL2_FOUND)
                message(WARNING "OpenBabel v2 not found.")
            endif()
        endif()
    endif()

    if(NOT (OPENBABEL2_FOUND OR OPENBABEL3_FOUND))
        if(NOT DOWNLOAD_DEPENDENCIES)
            message(FATAL_ERROR "OpenBabel is not available!")
        endif()

        message(STATUS "Will download and compile OpenBabel in place")
        set(OPENBABEL_LIB_FILE ${CMAKE_SOURCE_DIR}/external/openbabel-install/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openbabel${CMAKE_STATIC_LIBRARY_SUFFIX})
        ExternalProject_add(OpenBabel_external
            GIT_REPOSITORY  https://github.com/openbabel/openbabel.git
            GIT_TAG         openbabel-3-0-0
            GIT_SHALLOW     TRUE
            GIT_PROGRESS    TRUE
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/openbabel-src
            BINARY_DIR ${CMAKE_SOURCE_DIR}/external/openbabel-build
            CMAKE_ARGS -DBUILD_TESTING=OFF -DBUILD_MIXED=ON -DBUILD_SHARED=OFF
                       -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/openbabel-install
                       -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            BUILD_BYPRODUCTS ${OPENBABEL_LIB_FILE}
        )

        # Set openbabel variables manually
        set(OPENBABEL3_FOUND TRUE)
        set(OPENBABEL3_INCLUDE_DIR  ${CMAKE_SOURCE_DIR}/external/openbabel-install/include/openbabel3)
        set(OPENBABEL3_LIBRARIES    ${OPENBABEL_LIB_FILE})
    endif()
endif()

#--------------------
# Gromacs
#--------------------
if(WITH_GROMACS)
     # See if pathes are provided
     if(TRY_SYSTEM_GROMACS AND GROMACS_SRC_ROOT AND GROMACS_LIBRARIES)
         message(STATUS "Gromacs sources root set manually to ${GROMACS_SRC_ROOT}")
         message(STATUS "Gromacs libraries set manually to ${GROMACS_LIBRARIES}")
     else()
        if(NOT DOWNLOAD_DEPENDENCIES)
            message(FATAL_ERROR "Gromacs is not available!")
        endif()

        message(STATUS "Will download and compile Gromacs in place")

        FetchContent_Declare(Gromacs_external_fetch
            GIT_REPOSITORY  https://gitlab.com/gromacs/gromacs.git
            GIT_TAG         master #v2020.5
            GIT_SHALLOW     TRUE
            GIT_PROGRESS    TRUE
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/gromacs-src
        )
        FetchContent_GetProperties(Gromacs_external_fetch)
        if(NOT Gromacs_external_fetch_POPULATED)
          FetchContent_Populate(Gromacs_external_fetch)
        endif()

        set(GROMACS_LIB_FILE ${CMAKE_SOURCE_DIR}/external/gromacs-build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX})
        ExternalProject_add(Gromacs_external
            #GIT_REPOSITORY  https://gitlab.com/gromacs/gromacs.git
            #GIT_TAG         master #v2020.5
            #GIT_SHALLOW     TRUE
            #GIT_PROGRESS    TRUE
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/gromacs-src
            BINARY_DIR ${CMAKE_SOURCE_DIR}/external/gromacs-build
            CMAKE_ARGS  -DGMX_MPI=OFF -DGMX_GPU=OFF -DGMX_SIMD=none
                        -DGMX_FFT_LIBRARY=fftpack
                        -DBUILD_TESTING=OFF -DGMXAPI=OFF -DGMX_IMD=OFF
                        -DGMX_INSTALL_NBLIB_API=OFF -DGMX_OPENMP=OFF
                        -DGMX_THREAD_MPI=OFF
                        -DBUILD_SHARED_LIBS=OFF #-DGMX_USE_TNG=OFF
                        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            INSTALL_COMMAND ""
            BUILD_BYPRODUCTS ${GROMACS_LIB_FILE}
        )
        set(GROMACS_SRC_ROOT    ${CMAKE_SOURCE_DIR}/external/gromacs-src)
        set(GROMACS_LIBRARIES   ${GROMACS_LIB_FILE})
    endif()

    # Get actual gromacs version
    file(READ ${GROMACS_SRC_ROOT}/cmake/gmxVersionInfo.cmake gmxinfofile)
    string(REGEX MATCH "set\\(GMX_VERSION_MAJOR +([0-9]+)" match ${gmxinfofile})
    set(GROMACS_VERSION ${CMAKE_MATCH_1})
    # Now we have GMX_VERSION_MAJOR in GROMACS_VERSION!

    if(GROMACS_VERSION GREATER 2020)
        set(GROMACS_INCLUDE_DIRECTORIS
            ${CMAKE_SOURCE_DIR}/external/gromacs-src/src                  # Gromacs up to 2020.5
            ${CMAKE_SOURCE_DIR}/external/gromacs-src/api/legacy/include   # Gromacs 2021.x
            ${CMAKE_SOURCE_DIR}/external/gromacs-src/src/external         # Gromacs 2021.x
        )
    else()
        set(GROMACS_INCLUDE_DIRECTORIS
            ${CMAKE_SOURCE_DIR}/external/gromacs-src/src                  # Gromacs up to 2020.5
        )
    endif()

    # Configure include file with Gromacs version
    configure_file(${PROJECT_SOURCE_DIR}/src/core/gromacs_utils/gromacs_version_info.h.in
                   ${CMAKE_CURRENT_BINARY_DIR}/src/core/gromacs_utils/gromacs_version_info.h @ONLY)
endif()

#--------------------
# TNG_IO
#--------------------
if(WITH_TNG)
    if(NOT DOWNLOAD_DEPENDENCIES)
        message(FATAL_ERROR "tng_io is not available!")
    endif()

    message(STATUS "Will download and compile tng_io library in place")
    set(TNG_LIB_FILE ${CMAKE_SOURCE_DIR}/external/tng-install/lib/${CMAKE_STATIC_LIBRARY_PREFIX}tng_io${CMAKE_STATIC_LIBRARY_SUFFIX})
    ExternalProject_add(TNG_external
        GIT_REPOSITORY  https://gitlab.com/gromacs/tng.git
        GIT_TAG         v1.8.2
        GIT_SHALLOW     TRUE
        GIT_PROGRESS    TRUE
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/tng-src
        BINARY_DIR ${CMAKE_SOURCE_DIR}/external/tng-build
        CMAKE_ARGS  -DBUILD_SHARED_LIBS=OFF
                    -DTNG_BUILD_EXAMPLES=OFF
                    -DTNG_BUILD_TEST=OFF
                    -DTNG_BUILD_OWN_ZLIB=ON
                    -DBUILD_TESTING=OFF
                    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                    -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/tng-install
        BUILD_BYPRODUCTS ${TNG_LIB_FILE}
    )
    set(TNG_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/tng-install/include)
    set(TNG_LIBRARIES   ${TNG_LIB_FILE})
endif()
