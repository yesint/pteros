add_library(pteros_sim 
simulation.cpp 
gromacs_ff.cpp 
${PROJECT_SOURCE_DIR}/include/pteros/simulation/simulation.h
${PROJECT_SOURCE_DIR}/include/pteros/simulation/gromacs_ff.h)

# Link with pteros libs
link_directories(${PROJECT_SOURCE_DIR}/lib)
target_link_libraries(pteros_sim  pteros)

# Only need parent pteros include dir
include_directories(${PROJECT_SOURCE_DIR}/include)

set_target_properties(pteros_sim PROPERTIES
ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)



