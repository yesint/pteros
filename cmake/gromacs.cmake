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

if(WITH_GROMACS)
     # See if pathes are provided
     if(TRY_SYSTEM_GROMACS AND GROMACS_SOURCE_DIR AND GROMACS_BINARY_DIR)
         message(STATUS "Gromacs sources root set manually to ${GROMACS_SOURCE_DIR}")
         message(STATUS "Gromacs binary dir set manually to ${GROMACS_BINARY_DIR}")
     else()        
        message(STATUS "Will download and compile Gromacs")

        CPMAddPackage(
            NAME            GROMACS
            GIT_REPOSITORY  https://gitlab.com/gromacs/gromacs.git
            GIT_TAG         main
            DOWNLOAD_ONLY   ON
        )
        # Variable GROMACS_SOURCE_DIR and GROMACS_BINARY_DIR are set by CMP at this point

        # File to be set as a byproduct to trigger buld correctly
        set(GROMACS_LIB_FILES
            ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}muparser${CMAKE_STATIC_LIBRARY_SUFFIX}
        )

        # Get actual gromacs version
        file(READ ${GROMACS_SOURCE_DIR}/cmake/gmxVersionInfo.cmake gmxinfofile)
        string(REGEX MATCH "set\\(GMX_VERSION_MAJOR +([0-9]+)" match ${gmxinfofile})
        if(CMAKE_MATCH_1)
            set(GROMACS_VERSION ${CMAKE_MATCH_1})
        else()
            # Since Gromacs 2022.x version is stored in main CmakeLists file
            file(READ ${GROMACS_SOURCE_DIR}/CMakeLists.txt gmxinfofile)
            string(REGEX MATCH "Gromacs VERSION +([0-9]+)" match ${gmxinfofile})
            set(GROMACS_VERSION ${CMAKE_MATCH_1})
        endif()
        message(STATUS "Gromacs version used: ${GROMACS_VERSION}")

        # Building Gromacs without TNG is broken before 2023.1
        if(GROMACS_VERSION VERSION_GREATER_EQUAL 2023.1)
            set(GMX_USE_TNG_FLAG OFF)
        else()
            set(GMX_USE_TNG_FLAG ON)
        endif()

        # The only Gromacs subdirectories needed are fileio, mdtypes, topology
        # All the rest could be switched off
        #set(GMX_NEEDED fileio linearalgebra math mdtypes pbcutil simd topology utility)

        ExternalProject_add(Gromacs_external
            SOURCE_DIR ${GROMACS_SOURCE_DIR}
            BINARY_DIR ${GROMACS_BINARY_DIR}
            CMAKE_ARGS  -DGMX_MPI=OFF
                        -DGMX_GPU=OFF
                        -DGMX_SIMD=none
                        -DGMX_FFT_LIBRARY=fftpack
                        -DBUILD_TESTING=OFF
                        -DGMXAPI=OFF
                        -DGMX_IMD=OFF
                        -DGMX_INSTALL_NBLIB_API=OFF
                        -DGMX_OPENMP=OFF
                        -DGMX_THREAD_MPI=OFF
                        -DBUILD_SHARED_LIBS=OFF
                        -DGMX_USE_TNG=${GMX_USE_TNG_FLAG}
                        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            INSTALL_COMMAND ""
            BUILD_BYPRODUCTS ${GROMACS_LIB_FILES}
        )
        set(GROMACS_LIBRARIES ${GROMACS_LIB_FILES})
    endif()    

    if(GROMACS_VERSION VERSION_GREATER 2022) #2023 and up
        set(GROMACS_INCLUDE_DIRECTORIS
            ${GROMACS_SOURCE_DIR}/src                  # Gromacs up to 2020.5
            ${GROMACS_SOURCE_DIR}/src/gromacs/utility/include # Gromacs 2023.x
            ${GROMACS_SOURCE_DIR}/src/gromacs/math/include # Gromacs 2023.x
            ${GROMACS_SOURCE_DIR}/src/gromacs/topology/include # Gromacs 2023.x
            ${GROMACS_SOURCE_DIR}/api/legacy/include   # Gromacs 2021.x
            ${GROMACS_BINARY_DIR}/api/legacy/include   # Gromacs 2023.x
            ${GROMACS_SOURCE_DIR}/src/external         # Gromacs 2021.x
        )
        set(GROMACS_LIBRARIES
           ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX}
           ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}muparser${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
    elseif(GROMACS_VERSION VERSION_GREATER 2020) #2021-2022
        set(GROMACS_INCLUDE_DIRECTORIS
            ${GROMACS_SOURCE_DIR}/src                  # Gromacs up to 2020.5
            ${GROMACS_SOURCE_DIR}/api/legacy/include   # Gromacs 2021.x
            ${GROMACS_SOURCE_DIR}/src/external         # Gromacs 2021.x
        )
        set(GROMACS_LIBRARIES 
            ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX}
            ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}muparser${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
    else() #2020 and below
        set(GROMACS_INCLUDE_DIRECTORIS
            ${GROMACS_SOURCE_DIR}/src                  # Gromacs up to 2020.5
        )
        set(GROMACS_LIBRARIES 
            ${GROMACS_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gromacs${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
    endif()

    # Configure include file with Gromacs version
    configure_file(${PROJECT_SOURCE_DIR}/src/core/gromacs_version_info.h.in
                   ${CMAKE_BINARY_DIR}/src/core/gromacs_version_info.h @ONLY)

    #--------------------------------------------------------------------------------
    # Create a Gromacs interface library to provide headers and libs to other targets
    #--------------------------------------------------------------------------------
    add_library(gromacs_interface INTERFACE)
    target_include_directories(gromacs_interface INTERFACE
        ${GROMACS_INCLUDE_DIRECTORIS}
        ${CMAKE_BINARY_DIR}/src/core/  # For generated gromacs_version_info.h
        )
    target_link_libraries(gromacs_interface INTERFACE ${GROMACS_LIBRARIES})
    target_compile_definitions(gromacs_interface INTERFACE USE_GROMACS)
    if(TARGET Gromacs_external)
        add_dependencies(gromacs_interface Gromacs_external)
    endif()
endif()
