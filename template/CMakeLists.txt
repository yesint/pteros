###############################################
## This is general CMakeLists.txt file       ##
## which is needed to link your own program  ##
## with Pteros libraries                     ##
###############################################
cmake_minimum_required(VERSION 2.6)
project(pteros_template_program)

###############################################
## You have to set the variables below:      ##
###############################################

SET(PTEROS_INCLUDE_PATH "/home/yevhen/programs/pteros/include")
SET(PTEROS_LIB_PATH "/home/yevhen/programs/pteros/lib")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    ##################################################
    ## On Windows put the path to boost libs here:  ##
    ##################################################

    # Path to Boost:
    SET(BOOST_ROOT "c:/install/boost")
    # Path to Boost includes:
    SET(BOOST_INCLUDEDIR "C:/install/boost/include/boost-1_53")
    # Path to Boost libs:
    link_directories("C:/install/boost/lib")
ENDIF()


#############################################################################
## !!! Do not modify anything below unless you know CMake well enough !!!  ##
#############################################################################

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -msse2 -std=c++11")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")    
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_THREAD_USE_LIB -msse2 -std=c++11")
ENDIF()

# Strip libraries and executables to decrease size in Release builds
SET(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -s")
SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s")


SET(CMAKE_BUILD_TYPE Release)

find_package(Boost 1.48.0 REQUIRED system date_time)

include_directories(${PTEROS_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})
link_directories(${PTEROS_LIB_PATH})
link_libraries(pteros pteros_analysis ${Boost_LIBRARIES})

add_subdirectory(src)
