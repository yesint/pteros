###############################################
## This is CMakeLists.txt file               ##
## for Pteros compiled analysis plugin       ##
###############################################
cmake_minimum_required(VERSION 2.6)
project(pteros_analysis_plugin)

###############################################
## You have to set the variables below:      ##
###############################################

# Name of your plugin. This becomes the name of Python module!
SET(PLUGIN_NAME example_plugin)

# The list of source files for this plugin
SET(PLUGIN_FILES plugin_main.cpp)

# Path to your Pteros installation
SET(PTEROS_PATH /home/semen/programs/pteros/)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    ##################################################################
    ## On Windows put the full path to boost and python libs here:  ##
    ##################################################################

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
ENDIF()

#############################################################################
## !!! Do not modify anything below unless you know what are you doing !!! ##
#############################################################################

SET(PTEROS_INCLUDE_PATH "${PTEROS_PATH}/include")
SET(PTEROS_LIB_PATH "${PTEROS_PATH}/lib")

SET(CMAKE_INSTALL_PREFIX ${PTEROS_PATH})

add_definitions(-DPLUGIN_NAME=${PLUGIN_NAME})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(LIB_EXT ".pyd") # Force .pyd extension on Windows to ensure correct loading
ELSE()
    SET(LIB_EXT ".so") # Use default on Linux
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -msse2")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_THREAD_USE_LIB -msse2 -s -std=c++0x")
ENDIF()

# Strip libraries and executables to decrease size in Release builds
SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -s")
SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s")

SET(CMAKE_BUILD_TYPE Release)

# Will use shared Boost libraries
set(Boost_USE_STATIC_LIBS OFF)

find_package(Boost 1.48.0 REQUIRED regex system date_time python)
find_package(PythonInterp)
find_package(PythonLibs)

include_directories(${PTEROS_INCLUDE_PATH} ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
link_directories(${PTEROS_LIB_PATH})
link_libraries(pteros pteros_analysis ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})


add_library(plugin SHARED ${PLUGIN_FILES})

set_target_properties(plugin PROPERTIES
    PREFIX ""
    SUFFIX "${LIB_EXT}"
    OUTPUT_NAME "${PLUGIN_NAME}"
)

#Installation
install(TARGETS plugin
    LIBRARY DESTINATION python/pteros_analysis_plugins
    ARCHIVE DESTINATION python/pteros_analysis_plugins
)
