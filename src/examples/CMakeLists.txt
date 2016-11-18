link_libraries(pteros pteros_analysis)

add_executable(parallel_processing parallel_processing.cpp)
add_executable(benchmark benchmark.cpp)
add_executable(cg_lipids cg_lipids.cpp)
add_executable(chol_counter chol_counter.cpp bilayer.h bilayer.cpp)
add_executable(water_movement water_movement.cpp)

install(TARGETS parallel_processing benchmark cg_lipids chol_counter water_movement
    RUNTIME DESTINATION bin/examples
)
