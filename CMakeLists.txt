cmake_minimum_required(VERSION 3.1)

project(pteros)

OPTION(PYTHON_BINDINGS "Build python bindings and compiled analysis plugins as Python modules" ON)
OPTION(POWERSASA "Use POWERSASA code. This implies severe license restrictions! \n Read the comment in src/core/sasa/power_sasa.h carefully before enabling this option!" OFF)
OPTION(STANDALONE_PLUGINS "Compile analysis plugins as stand-alone executables" OFF)
OPTION(MAKE_PACKAGE "Generate package" OFF)
OPTION(USE_OPENMP "Use OpenMP parallelization if supported" ON)

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

if(USE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Find Python only if we need bindings
if(PYTHON_BINDINGS)
    find_package(PythonInterp)
    find_package(PythonLibs)    

    if(PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
        SET(WITH_PYTHON 1) # Compile python wrappers
        #For conditional compilation
        add_definitions(-DWRAP_PYTHON)

        # Set python includes globally
        include_directories(${PYTHON_INCLUDE_PATH})

        # Find boost.python
        find_package(Boost 1.48.0 REQUIRED python)
        # Save lib location separately from other boost libs
        set(BOOST_PYTHON_LIB ${Boost_LIBRARIES})

        #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION")
    else()
        message("Can't find Python interpreter of libraries!\nPython bindings and plugins are disabled!")
    endif()
endif()

# Find Boost
find_package(Boost 1.48.0 REQUIRED system date_time)

# For boost>1.50 and on windows also add filesystem to compile normally
if(${Boost_MINOR_VERSION} GREATER 50 AND ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    find_package(Boost 1.48.0 REQUIRED filesystem system)
endif()

# Add Boost includes and libs globally - they are needed everywhere
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# Add Eigen
find_package(Eigen3 3.2 QUIET)
if(NOT EIGEN3_INCLUDE_DIR)
  set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/thirdparty)
  install(DIRECTORY ${EIGEN3_INCLUDE_DIR}/Eigen DESTINATION include)
  message("Eigen3 library was not found in the system. Falling back to ${EIGEN3_INCLUDE_DIR}.")
endif()
include_directories(${EIGEN3_INCLUDE_DIR})

# Add fmt library support
#add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/fmt EXCLUDE_FROM_ALL)
#include_directories(${PROJECT_SOURCE_DIR}/thirdparty/fmt)
#for spdlog
include_directories(${PROJECT_SOURCE_DIR}/thirdparty)

# Directories for compiled files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Main Pteros include dir (for build time, not yet installed)
include_directories(include)

# Main Pteros lib dir (for build time, not yet installed)
link_directories("${CMAKE_BINARY_DIR}/lib")

# Determine extension of the shared libs depending on platform
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(SHARED_LIB_EXT ".pyd") # Force .pyd extension on Windows to ensure correct loading
ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(SHARED_LIB_EXT ".so") # Use .so on Linux
ENDIF()

# Set python install dir
set(PY_INST_DIR python)

#--------------
# CPack config
#--------------
IF(MAKE_PACKAGE)
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Pteros molecular modeling library")
    SET(CPACK_PACKAGE_VERSION_MAJOR "2")
    SET(CPACK_PACKAGE_VERSION_MINOR "0")
    SET(CPACK_PACKAGE_VERSION_PATCH "0")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        SET(CPACK_GENERATOR DEB)
        SET(CPACK_PACKAGE_CONTACT "Semen Yesylevskyy yesint4@yahoo.com")
        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "g++ (>=4.8), libboost-all-dev (>=1.50), python-dev (>=2.5), python-numpy")
        execute_process(COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(PY_INST_DIR ${PYTHON_SITE_PACKAGES})
    ENDIF()
    INCLUDE(CPack)
ENDIF()
#--------------

add_subdirectory(src/core)
add_subdirectory(src/core/sasa)
add_subdirectory(src/analysis)
add_subdirectory(src/tools)
add_subdirectory(src/python)
add_subdirectory(src/python/bindings)
add_subdirectory(src/python/compiled_analysis_plugins)
add_subdirectory(src/examples)
add_subdirectory(src/test)
#add_subdirectory(template)
#add_subdirectory(template_plugin)

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
