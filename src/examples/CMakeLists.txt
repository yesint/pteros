SET(PTEROS_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/include")
SET(PTEROS_LIB_PATH "${PROJECT_SOURCE_DIR}/lib")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -msse2") 
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC") 

SET(CMAKE_BUILD_TYPE Release)

include_directories(${PTEROS_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})
link_directories(${PTEROS_LIB_PATH})
link_libraries(pteros pteros_analysis pteros_sim xdrfile ${Boost_LIBRARIES})

add_executable(parallel_processing parallel_processing.cpp)
add_executable(benchmark benchmark.cpp)
add_executable(cg_lipids cg_lipids.cpp)
add_executable(chol_counter chol_counter.cpp)
add_executable(water_movement water_movement.cpp)
