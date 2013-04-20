set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples")

link_libraries(pteros pteros_analysis ${Boost_LIBRARIES})

add_executable(parallel_processing parallel_processing.cpp)
add_executable(benchmark benchmark.cpp)
add_executable(cg_lipids cg_lipids.cpp)
add_executable(chol_counter chol_counter.cpp)
add_executable(water_movement water_movement.cpp)
