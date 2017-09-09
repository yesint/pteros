if(POWERSASA)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
include_directories(${PROJECT_SOURCE_DIR}/src/core/sasa)

add_executable(sasatest
${PROJECT_SOURCE_DIR}/src/core/sasa/sasatest.cpp
${PROJECT_SOURCE_DIR}/src/core/sasa/power_diagram.h
${PROJECT_SOURCE_DIR}/src/core/sasa/power_sasa.h
)

endif()
