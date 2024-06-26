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

if(WITH_OPENBABEL)
    if(TRY_SYSTEM_OPENBABEL)
        find_package(OpenBabel3 3.0.0)
     endif()

    if(NOT OPENBABEL3_FOUND)
        message(STATUS "Will download and compile OpenBabel")
        CPMAddPackage(
            NAME            OPENBABEL
            GIT_REPOSITORY  https://github.com/openbabel/openbabel.git
            GIT_TAG         openbabel-3-0-0
            DOWNLOAD_ONLY   ON
        )
        # Variable OPENBABEL_SOURCE_DIR and OPENBABEL_BINARY_DIR are set by CMP at this point

        set(OPENBABEL_INSTALL_DIR ${CMAKE_BINARY_DIR}/openbabel-install)
        # File to be set as a byproduct to trigger buld correctly
        set(OPENBABEL_LIB_FILE ${OPENBABEL_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openbabel${CMAKE_STATIC_LIBRARY_SUFFIX})

        ExternalProject_add(OpenBabel_external
            SOURCE_DIR ${OPENBABEL_SOURCE_DIR}
            BINARY_DIR ${OPENBABEL_BINARY_DIR}
            CMAKE_ARGS -DBUILD_TESTING=OFF
                       -DBUILD_MIXED=ON
                       -DBUILD_SHARED=OFF
                       -DCMAKE_INSTALL_PREFIX=${OPENBABEL_INSTALL_DIR}
                       -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            BUILD_BYPRODUCTS ${OPENBABEL_LIB_FILE}
        )

        # Set openbabel variables
        set(OPENBABEL3_FOUND TRUE)
        set(OPENBABEL3_INCLUDE_DIR
            ${OPENBABEL_INSTALL_DIR}/include/openbabel3
            ${OPENBABEL_SOURCE_DIR}/include/openbabel3
        )
        set(OPENBABEL3_LIBRARIES ${OPENBABEL_LIB_FILE})
    endif()

    #--------------------------------------------------------------------------------
    # Create a Babel interface library to provide headers and libs to other targets
    #--------------------------------------------------------------------------------
    add_library(openbabel_interface INTERFACE)
    target_include_directories(openbabel_interface INTERFACE ${OPENBABEL3_INCLUDE_DIR})
    target_link_libraries(openbabel_interface INTERFACE ${OPENBABEL3_LIBRARIES})
    target_compile_definitions(openbabel_interface INTERFACE USE_OPENBABEL)
    if(TARGET OpenBabel_external)
        add_dependencies(openbabel_interface OpenBabel_external)
    endif()
endif()
