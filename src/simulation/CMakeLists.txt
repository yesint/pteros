add_library(pteros_sim SHARED
simulation.cpp 
gromacs_ff.cpp 
${PROJECT_SOURCE_DIR}/include/pteros/simulation/simulation.h
${PROJECT_SOURCE_DIR}/include/pteros/simulation/gromacs_ff.h)

# Only need parent pteros include dir
include_directories(${PROJECT_SOURCE_DIR}/include)

#Installation
install(TARGETS pteros_sim
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)



