cmake_minimum_required(VERSION 3.1)

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


###############################################
## Do not modify below!                      ##
###############################################

set_source_files_properties(${PLUGIN_FILES} PROPERTIES COMPILE_DEFINITIONS PLUGIN_NAME=${PLUGIN_NAME})

SET(PTEROS_INCLUDE_PATH "${PTEROS_PATH}/include")
SET(PTEROS_LIB_PATH "${PTEROS_PATH}/lib")
SET(CMAKE_INSTALL_PREFIX ${PTEROS_PATH})

include_directories(${PTEROS_INCLUDE_PATH})
link_directories(${PTEROS_LIB_PATH})

# Gloabal compiler and linker flags for Linux
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    # Strip libraries and executables to decrease size in Release builds
    SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -s")
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s")
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

# Find Python
message("Looking for Python...")
find_package(pybind11)
# Set python includes globally (pybind11 fails to find Python.h without this)
include_directories(${PYTHON_INCLUDE_PATH})

# Find Eigen
message("Looking for Eigen...")
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

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


# Make python extension module
pybind11_add_module(plugin ${PLUGIN_FILES})

target_link_libraries(plugin PRIVATE pteros pteros_analysis)

set_target_properties(plugin PROPERTIES
    OUTPUT_NAME ${PLUGIN_NAME}
)

#Installation
install(TARGETS plugin
    LIBRARY DESTINATION python/pteros_analysis_plugins
    ARCHIVE DESTINATION python/pteros_analysis_plugins
)
