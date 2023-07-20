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

# Path to CMake modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# OpenMP
if(WITH_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
    if(NOT OpenMP_CXX_FOUND)
        message(WARNING "OpenMP is not available.")
    endif()
endif()

#======================================================
# Dependencies which can use CPM for installing
#======================================================

#--------------------
# Eigen
#--------------------
CPMAddPackage(
    NAME              Eigen
    GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
    GIT_TAG           ${EIGEN_VERSION}
    OPTIONS
        "EIGEN_BUILD_DOC OFF"
        "BUILD_TESTING OFF"
        "EIGEN_BUILD_TESTING OFF"
        "EIGEN_BUILD_PKGCONFIG OFF"
)

#--------------------
# fmt
#--------------------
CPMAddPackage(
    NAME                fmt
    GITHUB_REPOSITORY   fmtlib/fmt
    GIT_TAG             ${FMT_VERSION}
    OPTIONS
        "FMT_INSTALL ON"
        "FMT_DOC OFF"
        "FMT_TEST OFF"
)

#--------------------
# spdlog
#--------------------
CPMAddPackage(
    NAME                spdlog
    GITHUB_REPOSITORY   gabime/spdlog
    GIT_TAG             v${SPDLOG_VERSION}
    OPTIONS
        "SPDLOG_INSTALL ON"
        "SPDLOG_BUILD_TESTS OFF"
        "SPDLOG_BUILD_EXAMPLE OFF"
        "SPDLOG_FMT_EXTERNAL ON"
)

#--------------------
# Pybind11
#--------------------
if(WITH_PYTHON)    
    # Force Python3
    set(PYBIND11_PYTHON_VERSION 3 CACHE INTERNAL "" FORCE)

    CPMAddPackage(
        NAME                pybind11
        GITHUB_REPOSITORY   pybind/pybind11
        GIT_TAG             ${PYBIND11_VERSION}
    )

    # Manage python modules installation paths (yes, this is a mess)
    if(NOT DEFINED PY_INST_DIR)
        if(DEFINED ENV{VIRTUAL_ENV})
            # This gives relative path to the python packages
            # It provides correct path for venv, but in global
            # context it returns global dist-packages, so shouldn't be used
            execute_process(
                COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
                from distutils import sysconfig as sc
                print(sc.get_python_lib(prefix='', plat_specific=True))"
                OUTPUT_VARIABLE PY_INST_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            # This gives a relative path, which will be combined with CMAKE_INSTALL_PREFIX
        else()
            # We are not in venv, so install to site-packages, not to dist-packages
            execute_process(
                COMMAND "${PYTHON_EXECUTABLE}" -m site --user-site
                OUTPUT_VARIABLE PY_INST_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            # This gives correct absolute path directly
        endif()
    endif()
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
