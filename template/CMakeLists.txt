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

SET(PTEROS_INCLUDE_PATH "/home/semen/programs/pteros/include")
SET(PTEROS_LIB_PATH "/home/semen/programs/pteros/lib")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    ##################################################
    ## On Windows put the path to boost libs here:  ##
    ##################################################

    SET(BOOST_ROOT "c:/install/boost")
ENDIF()


#############################################################################
## !!! Do not modify anything below unless you know CMake well enough !!!  ##
#############################################################################

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -msse2")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_THREAD_USE_LIB -msse2")
ENDIF()

SET(CMAKE_BUILD_TYPE Release)

find_package(Boost 1.42.0 REQUIRED thread)

include_directories(${PTEROS_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})
link_directories(${PTEROS_LIB_PATH})
link_libraries(pteros pteros_analysis ${Boost_LIBRARIES})

add_subdirectory(src)
