###############################################
## This is CMakeLists.txt file               ##
## for Pteros compiled analysis plugin       ##
###############################################
cmake_minimum_required(VERSION 2.6)
project(example_analysis_plugin)

###############################################
## You have to set the variables below:      ##
###############################################

SET(PTEROS_PATH "/home/semen/programs/pteros/")
SET(PLUGIN_NAME "example_plugin")

# If you put Pteros include and lib directories elsewhere specify them here explicitly
SET(PTEROS_INCLUDE_PATH "${PTEROS_PATH}/include")
SET(PTEROS_LIB_PATH "${PTEROS_PATH}/lib")
SET(CMAKE_INSTALL_PREFIX ${PTEROS_PATH})

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    #############################################################
    ## On Windows put the path to boost and python libs here:  ##
    #############################################################

    # Path to Boost:
    SET(BOOST_ROOT "c:/install/boost")
    SET(BOOST_INCLUDEDIR "C:/install/boost/include/boost-1_53")
    link_directories("C:/install/boost/lib")

    # Path to numpy:
    include_directories("C:/Python27/Lib/site-packages/numpy/core/include")
    link_directories("C:/Python27/Lib/site-packages/numpy/core/lib")
ENDIF()


#############################################################################
## !!! Do not modify anything below unless you know CMake well enough !!!  ##
#############################################################################

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(LIB_EXT ".pyd") # Force .pyd extension on Windows to ensure correct loading
ELSE()
    SET(LIB_EXT ".so") # Use default on Linux
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -msse2 -s -std=c++0x")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -s")
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_THREAD_USE_LIB -msse2 -s -std=c++0x")
ENDIF()

SET(CMAKE_BUILD_TYPE Release)

# Will use shared Boost libraries
set(Boost_USE_STATIC_LIBS OFF)

find_package(Boost 1.46.0 REQUIRED thread regex python)
find_package(PythonInterp)
find_package(PythonLibs)

include_directories(${PTEROS_INCLUDE_PATH} ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
link_directories(${PTEROS_LIB_PATH})
link_libraries(pteros pteros_analysis ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})


add_library(plugin SHARED 
    "${PLUGIN_NAME}.cpp"
)

set_target_properties(plugin
PROPERTIES
PREFIX ""
SUFFIX "${LIB_EXT}"
OUTPUT_NAME "${PLUGIN_NAME}"
)

#Installation
install(TARGETS plugin
    LIBRARY DESTINATION python/pteros_analysis_plugins
    ARCHIVE DESTINATION python/pteros_analysis_plugins
)
