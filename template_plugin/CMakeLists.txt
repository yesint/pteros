cmake_minimum_required(VERSION 3.1)

project(pteros_analysis_plugin)

find_package(pybind11)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
find_package(PythonInterp  3.5)
find_package(PythonLibs 3.5)




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

# Find Python
message("Looking for Python...")
find_package(PythonInterp)
find_package(PythonLibs)
if(PYTHONINTERP_FOUND AND PYTHONLIBS_FOUND)
    # Configure pybind11
    message("Configuring Pybind11...")
    find_package(pybind11)
    # Set python includes globally (pybind11 fails to find Python.h without this)
    include_directories(${PYTHON_INCLUDE_PATH})
else()
    message("Can't find Python interpreter or libraries!")
endif()

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
