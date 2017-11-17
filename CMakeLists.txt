cmake_minimum_required(VERSION 3.1)

project(pteros)

OPTION(WITH_PYTHON "Build python bindings and compile analysis plugins as Python modules" ON)
OPTION(WITH_OPENMP "Use OpenMP parallelization if supported" ON)
OPTION(WITH_POWERSASA "Use POWERSASA code. This implies severe license restrictions! \n Read the comment in thirdparty/sasa/LICENSE carefully before enabling this option!" OFF)
OPTION(STANDALONE_PLUGINS "Compile analysis plugins as stand-alone executables" OFF)
OPTION(MAKE_PACKAGE "Generate package" OFF)
OPTION(INSTALL_INTERNAL_SPDLOG "Install internal spdlog headers in the case if there is no system-wide spdlog" ON)
OPTION(MAKE_EXAMPLES "Compile examples and plugin templates" OFF)


# Use shared Boost libraries
set(Boost_USE_STATIC_LIBS OFF)

# Gloabal compiler and linker flags for Linux
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    # Strip libraries and executables to decrease size in Release builds
    SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -s")
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s")
ENDIF()

# Windows-related stuff
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_THREAD_USE_LIB -s")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -s")

    ##########################################################
    # Cnange the paths here if you are compiling on Windows! #
    ##########################################################

    # Path to Boost:
    SET(BOOST_ROOT "c:/install/boost")
    # Path to Boost includes:
    SET(BOOST_INCLUDEDIR "C:/install/boost/include/boost-1_53")
    # Path to Boost libs:
    link_directories("C:/install/boost/lib")

    # Path to numpy includes:
    include_directories("C:/Python27/Lib/site-packages/numpy/core/include")
    # Path to numpy libs:
    link_directories("C:/Python27/Lib/site-packages/numpy/core/lib")

    ##########################################################
    # End of user-modifieble section!                        #
    ##########################################################
ENDIF()

#Sets optmization for Eigen in Release build
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse2")

# Enable STL debugging in debug mode
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")

# Enable C++11 and C99 everywhere
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

# Set specific options for GCC if used
if(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
endif()

# Find OpenMP if needed
if(WITH_OPENMP)
    message("Looking for OpenMP...")
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Find Python only if we need bindings
if(WITH_PYTHON)
    #For conditional compilation
    add_definitions(-DWRAP_PYTHON)
    # Configure pybind11
    message("Configuring Pybind11...")
    find_package(pybind11 REQUIRED)
    # Set python install dir
    set(PY_INST_DIR python)
endif()

# Find Boost
message("Looking for Boost...")
find_package(Boost 1.48.0 REQUIRED system date_time)
# For boost>1.50 and on windows also add filesystem to compile normally
if(${Boost_MINOR_VERSION} GREATER 50 AND ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    find_package(Boost 1.48.0 REQUIRED filesystem system)
endif()

# Add Boost includes and libs globally - they are needed everywhere
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# Add Eigen
message("Looking for Eigen...")
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


# Main Pteros include dir (for build time, not yet installed)
include_directories(include)


#----------------------------
# Build thirdparty libraries
#----------------------------

message("Looking for spdlog...")
find_package(spdlog CONFIG QUIET)
if(NOT ${spdlog_FOUND})
    message("spdlog not found! Using internal one.")
    include_directories(${PROJECT_SOURCE_DIR}/thirdparty)
    # Install spdlog headers
    if(${INSTALL_INTERNAL_SPDLOG})
        install(DIRECTORY ${PROJECT_SOURCE_DIR}/thirdparty/spdlog DESTINATION include)
    endif()
endif()

#include_directories for spdlog

# Build separate tng_io library (will be linked later)
add_subdirectory(thirdparty/tng_io)
include_directories(thirdparty/tng_io/include)

# DSSP library
add_subdirectory(thirdparty/dssp)
include_directories(thirdparty/dssp)

# VMD molfile plugins
add_subdirectory(thirdparty/molfile_plugins)
include_directories(thirdparty/molfile_plugins)

# xdrfile library
add_subdirectory(thirdparty/xdrfile)
include_directories(thirdparty/xdrfile)

#---------------------------------------------
# Compiling pteros itself
#---------------------------------------------

add_subdirectory(src/core)
add_subdirectory(src/analysis)
add_subdirectory(src/tools)
add_subdirectory(src/python)
add_subdirectory(src/python/compiled_analysis_plugins)
add_subdirectory(src/test)
add_subdirectory(src/extras)

IF(MAKE_EXAMPLES)
    add_subdirectory(src/examples)
    add_subdirectory(template_plugin)
ENDIF()

# Install headers
install(DIRECTORY include DESTINATION .)

#--------------
# CPack config
#--------------
IF(MAKE_PACKAGE)
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Pteros molecular modeling library")
    SET(CPACK_PACKAGE_VERSION_MAJOR "2")
    SET(CPACK_PACKAGE_VERSION_MINOR "5")
    SET(CPACK_PACKAGE_VERSION_PATCH "0")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        SET(CPACK_GENERATOR DEB)
        SET(CPACK_PACKAGE_CONTACT "Semen Yesylevskyy yesint4@yahoo.com")
        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "g++ (>=4.8), libboost-all-dev (>=1.50), python-dev (>=2.5), python-numpy, libeigen3-dev")
        execute_process(COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(PY_INST_DIR ${PYTHON_SITE_PACKAGES})
    ENDIF()
    INCLUDE(CPack)
ENDIF()


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
