link_libraries(pteros pteros_analysis)

add_executable(benchmark benchmark.cpp)
add_executable(cg_lipids cg_lipids.cpp)
add_executable(chol_counter chol_counter.cpp bilayer.h bilayer.cpp)
add_executable(water_movement water_movement.cpp)

install(TARGETS benchmark cg_lipids chol_counter water_movement
    RUNTIME DESTINATION share/pteros/examples
)
