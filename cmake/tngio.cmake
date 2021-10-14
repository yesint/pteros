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


if(WITH_TNG)
    if(NOT DOWNLOAD_DEPENDENCIES)
        message(FATAL_ERROR "tng_io is not available!")
    endif()

    message(STATUS "Will download and compile tng_io library in place")

    FetchContent_Declare(TNG_external_fetch
        GIT_REPOSITORY  https://gitlab.com/gromacs/tng.git
        GIT_TAG         v1.8.2
        GIT_SHALLOW     TRUE
        GIT_PROGRESS    TRUE
    )
    FetchContent_GetProperties(TNG_external_fetch)
    if(NOT TNG_external_fetch_POPULATED)
      FetchContent_Populate(TNG_external_fetch)
    endif()

    FetchContent_GetProperties(
        TNG_external_fetch
        SOURCE_DIR TNG_SOURCE_DIR
        BINARY_DIR TNG_BINARY_DIR
    )
    set(TNG_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/tng-install)
    set(TNG_LIB_FILE ${TNG_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}tng_io${CMAKE_STATIC_LIBRARY_SUFFIX})

    ExternalProject_add(TNG_external
        SOURCE_DIR ${TNG_SOURCE_DIR}
        BINARY_DIR ${TNG_BINARY_DIR}
        CMAKE_ARGS  -DBUILD_SHARED_LIBS=OFF
                    -DTNG_BUILD_EXAMPLES=OFF
                    -DTNG_BUILD_TEST=OFF
                    -DTNG_BUILD_OWN_ZLIB=ON
                    -DBUILD_TESTING=OFF
                    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                    -DCMAKE_INSTALL_PREFIX=${TNG_INSTALL_DIR}
                    -DCMAKE_INSTALL_LIBDIR=lib
        BUILD_BYPRODUCTS ${TNG_LIB_FILE}
    )
    set(TNG_INCLUDE_DIR ${TNG_INSTALL_DIR}/include ${TNG_SOURCE_DIR}/include)
    set(TNG_LIBRARIES   ${TNG_LIB_FILE})
endif()
