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


if(WITH_OPENBABEL)
    if(TRY_SYSTEM_OPENBABEL)        
        # Try to find OpenBabel 3
        find_package(OpenBabel3 3.0.0)
        if(NOT OPENBABEL3_FOUND)
            # Try to find OpenBabel 2
            message(STATUS "OpenBabel v3 not found, searching for v2...")
            find_package(OpenBabel2 2.4.9)
            if(NOT OPENBABEL2_FOUND)
                message(STATUS "OpenBabel v2 not found.")
            endif()
        endif()
    endif()

    if(NOT (OPENBABEL2_FOUND OR OPENBABEL3_FOUND))
        if(NOT DOWNLOAD_DEPENDENCIES)
            message(FATAL_ERROR "OpenBabel is not available!")
        endif()

        FetchContent_Declare(OpenBabel_external_fetch
            GIT_REPOSITORY  https://github.com/openbabel/openbabel.git
            GIT_TAG         openbabel-3-0-0
            GIT_SHALLOW     TRUE
            GIT_PROGRESS    TRUE
        )
        FetchContent_GetProperties(OpenBabel_external_fetch)
        if(NOT OpenBabel_external_fetch_POPULATED)
          FetchContent_Populate(OpenBabel_external_fetch)
        endif()

        FetchContent_GetProperties(
            OpenBabel_external_fetch
            SOURCE_DIR OPENBABEL_SOURCE_DIR
            BINARY_DIR OPENBABEL_BINARY_DIR
        )
        set(OPENBABEL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/openbabel-install)
        set(OPENBABEL_LIB_FILE ${OPENBABEL_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openbabel${CMAKE_STATIC_LIBRARY_SUFFIX})

        message(STATUS "Will download and compile OpenBabel in place")
        ExternalProject_add(OpenBabel_external
            SOURCE_DIR ${OPENBABEL_SOURCE_DIR}
            BINARY_DIR ${OPENBABEL_BINARY_DIR}
            CMAKE_ARGS -DBUILD_TESTING=OFF -DBUILD_MIXED=ON -DBUILD_SHARED=OFF
                       -DCMAKE_INSTALL_PREFIX=${OPENBABEL_INSTALL_DIR}
                       -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            BUILD_BYPRODUCTS ${OPENBABEL_LIB_FILE}
        )

        # Set openbabel variables manually
        set(OPENBABEL3_FOUND TRUE)
        set(OPENBABEL3_INCLUDE_DIR
            ${OPENBABEL_INSTALL_DIR}/include/openbabel3
            ${OPENBABEL_SOURCE_DIR}/include/openbabel3
        )
        set(OPENBABEL3_LIBRARIES ${OPENBABEL_LIB_FILE})
    endif()
endif()
