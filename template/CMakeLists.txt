###############################################
## This is general CMakeLists.txt file       ##
## which is needed to link your own program  ##
## with Pteros libraries                     ##
###############################################
cmake_minimum_required(VERSION 3.1)
project(pteros_template_program)

###############################################
## You have to set the variables below:      ##
###############################################

SET(PTEROS_INCLUDE_PATH "/home/semen/programs/pteros/include")
SET(PTEROS_LIB_PATH "/home/semen/programs/pteros/lib")


#############################################################################
## !!! Do not modify anything below unless you know CMake well enough !!!  ##
#############################################################################

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
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

# Enable C++11 and C99 everywhere
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

# Set specific options for GCC if used
if(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
endif()

# Find Eigen
message("Looking for Eigen...")
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find Boost
message("Looking for Boost...")
find_package(Boost 1.48.0 REQUIRED system date_time)
# For boost>1.50 and on windows also add filesystem to compile normally
if(${Boost_MINOR_VERSION} GREATER 50 AND ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    find_package(Boost 1.48.0 REQUIRED filesystem)
endif()
# Add Boost includes and libs globally - they are needed everywhere
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

message("Looking for Gromacs...")
find_package(gromacs QUIET)
if(NOT ${gromacs_FOUND})
    # Try /usr/local
    message("   Trying /usr/local...")
    set(gromacs_DIR "/usr/local/gromacs/share/cmake/gromacs")
    find_package(gromacs QUIET)
endif()

if(${gromacs_FOUND})
    message("-- Found Gromacs ${gromacs_VERSION} at ${GROMACS_INCLUDE_DIRS}")
    add_definitions(-DUSE_GROMACS)
    include_directories(${GROMACS_INCLUDE_DIRS})
    link_libraries(${GROMACS_LIBRARIES})
else()
    message("-- Gromacs not found!")
endif()

link_libraries(pteros pteros_analysis ${Boost_LIBRARIES})

add_subdirectory(src)
