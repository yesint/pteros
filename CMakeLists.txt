cmake_minimum_required(VERSION 3.1)

project(pteros)

OPTION(WITH_PYTHON "Build python bindings and compile analysis plugins as Python modules" ON)
OPTION(WITH_OPENMP "Use OpenMP parallelization if supported" ON)
OPTION(WITH_OPENBABEL "Use OpenBabel if supported" ON)
OPTION(WITH_GROMACS "Use Gromacs if supported. Required to read tpr files." ON)
OPTION(WITH_POWERSASA "Use POWERSASA code. This implies license restrictions described here: thirdparty/sasa/LICENSE" ON)
OPTION(MAKE_STANDALONE_PLUGINS "Compile analysis plugins as stand-alone executables" OFF)
OPTION(MAKE_PACKAGE "Generate package" OFF)
OPTION(INSTALL_INTERNAL_SPDLOG "Install internal spdlog headers in the case if there is no system-wide spdlog" ON)
OPTION(MAKE_EXAMPLES "Compile examples and plugin templates" OFF)
OPTION(MAKE_TEST "Compile tests" OFF)

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

##############################
# Unconditional dependencies:
##############################

# Boost
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost 1.50 REQUIRED COMPONENTS system)

# Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

##############################
# Conditional dependencies:
##############################

# OpenMP only in release mode
if(WITH_OPENMP AND CMAKE_BUILD_TYPE STREQUAL "Release")
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

# Python
if(WITH_PYTHON)    
    # Configure pybind11
    find_package(pybind11 REQUIRED)
    # Set python install dir
    set(PY_INST_DIR python)
endif()

# OpenBabel
if(WITH_OPENBABEL)
    set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
    find_package(OpenBabel2 REQUIRED)
endif()

# Gromacs
if(WITH_GROMACS)
    find_package(gromacs QUIET)
    if(NOT ${gromacs_FOUND})
        # Try /usr/local
        message(STATUS "Looking for Gromacs in /usr/local...")
        set(gromacs_DIR "/usr/local/gromacs/share/cmake/gromacs")
        find_package(gromacs REQUIRED)
    endif()

    if(${gromacs_FOUND})
        message("Found Gromacs ${gromacs_VERSION} at ${GROMACS_INCLUDE_DIRS}")
    endif()
endif()


#--------------
# CPack config
#--------------
IF(MAKE_PACKAGE)
    message("Configuring package...")
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Pteros molecular modeling library")
    SET(CPACK_PACKAGE_VERSION_MAJOR "2")
    SET(CPACK_PACKAGE_VERSION_MINOR "7")
    SET(CPACK_PACKAGE_VERSION_PATCH "0")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        SET(CPACK_GENERATOR DEB)
        SET(CPACK_PACKAGE_CONTACT "Semen Yesylevskyy yesint4@yahoo.com")
        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "g++ (>=4.8), libboost-all-dev (>=1.50), python-dev (>=2.7), python-numpy, libeigen3-dev")
        set(PY_INST_DIR ${PYTHON_SITE_PACKAGES})
    ENDIF()
    INCLUDE(CPack)
ENDIF()

#----------------------------
# Build thirdparty libraries
#----------------------------

# spdlog
find_package(spdlog CONFIG QUIET)
if(NOT ${spdlog_FOUND})
    message(STATUS "spdlog not found! Using internal one.")
    add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/spdlog)
    # Install spdlog headers
    if(${INSTALL_INTERNAL_SPDLOG})
        install(DIRECTORY ${PROJECT_SOURCE_DIR}/thirdparty/spdlog/include/spdlog
                DESTINATION include)
    endif()
endif()

# tng_io library
add_subdirectory(thirdparty/tng_io)

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

# Install headers
install(DIRECTORY include DESTINATION .)



IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    message("------------------------------")
    message("IMPORTANT!")
    message("In order to be able to link with Pteros and to use Pteros Python modules")
    message("you have to add the following lines to your .bashrc file:")
    message("export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib\"")
    message("export PYTHONPATH=\"\$PYTHONPATH:${CMAKE_INSTALL_PREFIX}/python\"")
    message("export PATH=\"\$PATH:${CMAKE_INSTALL_PREFIX}/bin\"")
    message("------------------------------")
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    message("------------------------------")
    message("IMPORTANT!")
    message("In order to be able to link with Pteros and to use Pteros Python modules")
    message("you have to add the following directories to your PATH:")
    message("${CMAKE_INSTALL_PREFIX}/lib")
    message("${CMAKE_INSTALL_PREFIX}/bin")
    message("${CMAKE_INSTALL_PREFIX}/python")
    message("------------------------------")
ENDIF()
