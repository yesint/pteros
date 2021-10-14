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


if(WITH_GROMACS)
     # See if pathes are provided
     if(TRY_SYSTEM_GROMACS AND GROMACS_SOURCE_DIR AND GROMACS_BINARY_DIR)
         message(STATUS "Gromacs sources root set manually to ${GROMACS_SOURCE_DIR}")
         message(STATUS "Gromacs binary dir set manually to ${GROMACS_BINARY_DIR}")
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
        )
        FetchContent_GetProperties(Gromacs_external_fetch)
        if(NOT Gromacs_external_fetch_POPULATED)
          FetchContent_Populate(Gromacs_external_fetch)
        endif()

        FetchContent_GetProperties(
            Gromacs_external_fetch
            SOURCE_DIR GROMACS_SOURCE_DIR
            BINARY_DIR GROMACS_BINARY_DIR
        )
        set(GROMACS_LIB_FILE ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX})

        ExternalProject_add(Gromacs_external
            SOURCE_DIR ${GROMACS_SOURCE_DIR}
            BINARY_DIR ${GROMACS_BINARY_DIR}
            CMAKE_ARGS  -DGMX_MPI=OFF -DGMX_GPU=OFF -DGMX_SIMD=none
                        -DGMX_FFT_LIBRARY=fftpack
                        -DBUILD_TESTING=OFF -DGMXAPI=OFF -DGMX_IMD=OFF
                        -DGMX_INSTALL_NBLIB_API=OFF -DGMX_OPENMP=OFF
                        -DGMX_THREAD_MPI=OFF
                        -DBUILD_SHARED_LIBS=OFF #-DGMX_USE_TNG=OFF
                        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            INSTALL_COMMAND ""
            BUILD_BYPRODUCTS ${GROMACS_LIB_FILE}
        )
        set(GROMACS_LIBRARIES ${GROMACS_LIB_FILE})
    endif()

    # Get actual gromacs version
    file(READ ${GROMACS_SOURCE_DIR}/cmake/gmxVersionInfo.cmake gmxinfofile)
    string(REGEX MATCH "set\\(GMX_VERSION_MAJOR +([0-9]+)" match ${gmxinfofile})
    set(GROMACS_VERSION ${CMAKE_MATCH_1})
    # Now we have GMX_VERSION_MAJOR in GROMACS_VERSION!

    if(GROMACS_VERSION GREATER 2020)
        set(GROMACS_INCLUDE_DIRECTORIS
            ${GROMACS_SOURCE_DIR}/src                  # Gromacs up to 2020.5
            ${GROMACS_SOURCE_DIR}/api/legacy/include   # Gromacs 2021.x
            ${GROMACS_SOURCE_DIR}/src/external         # Gromacs 2021.x
        )
	set(GROMACS_LIBRARIES 
  	   ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX}
	   ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}muparser${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
    else()
        set(GROMACS_INCLUDE_DIRECTORIS
            ${GROMACS_SOURCE_DIR}/src                  # Gromacs up to 2020.5
        )
	set(GROMACS_LIBRARIES 
           ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
    endif()

    # Configure include file with Gromacs version
    configure_file(${PROJECT_SOURCE_DIR}/src/core/gromacs_utils/gromacs_version_info.h.in
                   ${CMAKE_BINARY_DIR}/src/core/gromacs_utils/gromacs_version_info.h @ONLY)
endif()