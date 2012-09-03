cmake_minimum_required(VERSION 2.6)
project(pteros)

OPTION(DEBUG_MODE "Compile in debug mode")
OPTION(PYTHON_BINDINGS "Build python bindings" ON)
OPTION(CHAISCRIPT_BINDINGS "Build chaiscript bindings")
OPTION(VMD_PDB "Use VMD molfile plugin for reading PDB files (if OFF - internal reader is used)")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=gnu++0x")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC ")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x -DBOOST_THREAD_USE_LIB")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_THREAD_USE_LIB")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")   
    SET(PYTHON_BINDINGS OFF)
    message("Python bindings are disabled on Windows")
    SET(BOOST_ROOT "c:/install/boost")
ENDIF()

if(DEBUG_MODE)
    message("Debug build")
    SET(CMAKE_BUILD_TYPE Debug)
else()
    message("Release build")
    SET(CMAKE_BUILD_TYPE Release)
    #Sets optmization for Eigen
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2") 
endif()

if(PYTHON_BINDINGS)
    find_package(PythonInterp)
    find_package(PythonLibs)
    find_package(Boost 1.46.0 REQUIRED python)
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

find_package(Boost 1.46.0)
if(Boost_MINOR_VERSION LESS 49)
    find_package(Boost REQUIRED thread regex)
else()
    find_package(Boost REQUIRED thread regex system)
endif()

include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES}) 

include_directories(include)
link_libraries( xdrfile )

add_subdirectory(src/xdrfile)
add_subdirectory(src/core)
add_subdirectory(src/test)
add_subdirectory(src/simulation)
add_subdirectory(src/analysis)
add_subdirectory(src/python)
add_subdirectory(src/examples)
add_subdirectory(src/molfile_plugins)

if(CHAISCRIPT_BINDINGS)
    add_subdirectory(src/chaiscript)
endif()

