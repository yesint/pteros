cmake_minimum_required(VERSION 2.6)
project(pteros)

#SET(CMAKE_BUILD_TYPE Release)
OPTION(PYTHON_BINDINGS "Build python bindings" ON)
OPTION(VMD_PDB "Use VMD molfile plugin for reading PDB files (if OFF - internal reader is used)" ON)
OPTION(POWERSASA "Use POWERSASA code. This implies severe license restrictions! \n Read the comment in src/core/sasa/power_sasa.h carefully before enabling this option!" OFF)
OPTION(NO_CPP11 "Old compiler is used. Disable code with all C++11 features" OFF)

# Will use shared Boost libraries
set(Boost_USE_STATIC_LIBS OFF)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC ")
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_THREAD_USE_LIB -s")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -s")

    ##########################################################
    # Cnange the paths here if you are compiling on Windows! #
    ##########################################################

    # Path to Boost:
    SET(BOOST_ROOT "c:/install/boost")
    SET(BOOST_INCLUDEDIR "C:/install/boost/include/boost-1_53")
    link_directories("C:/install/boost/lib")

    # Path to numpy:
    include_directories("C:/Python27/Lib/site-packages/numpy/core/include")
    link_directories("C:/Python27/Lib/site-packages/numpy/core/lib")
ENDIF()

if(CMAKE_BUILD_TYPE STREQUAL Release)
    message("Release build")
    #Sets optmization for Eigen
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
else()
    message("Debug build")
endif()

if(PYTHON_BINDINGS)
    find_package(PythonInterp)
    find_package(PythonLibs)
    find_package(Boost 1.48.0 REQUIRED python)
    if(Boost_FOUND AND PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
        SET(WRAP_PYTHON 1) # Compile python wrappers
        #For conditional compile
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWRAP_PYTHON")

        include_directories(${PYTHON_INCLUDE_PATH})
        link_libraries(${PYTHON_LIBRARIES})
        # Latest test for boost was for boost::python module, so include it
        include_directories(${Boost_INCLUDE_DIRS})

        message("Building with python bindings")
    endif()
endif()

find_package(Boost 1.48.0)
if(Boost_MINOR_VERSION LESS 49)
    find_package(Boost REQUIRED thread regex date_time)
else()
    find_package(Boost REQUIRED thread regex system date_time)
endif()

# In debug build we also need boost.date_time to link correctly. I have no idea why...
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    find_package(Boost REQUIRED date_time)
endif()

message("${Boost_INCLUDE_DIRS}")
message( "${Boost_LIBRARIES}" )

include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

include_directories(include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
link_directories("${CMAKE_BINARY_DIR}/lib")


add_subdirectory(src/core)
add_subdirectory(src/test)
add_subdirectory(src/analysis)
add_subdirectory(src/python)
add_subdirectory(src/python/bindings)
add_subdirectory(src/python/compiled_analysis_plugins)
add_subdirectory(src/examples)
add_subdirectory(src/core/sasa)

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
