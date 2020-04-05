add_library(pteros_gromacs_utils STATIC
    gromacs_utils.cpp
    gromacs_utils.h
)

target_link_libraries(pteros_gromacs_utils
    PUBLIC
        Eigen3::Eigen
        pteros_interface
)

target_include_directories(pteros_gromacs_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(WITH_GROMACS)
    target_compile_definitions(pteros_gromacs_utils PRIVATE USE_GROMACS)
    target_include_directories(pteros_gromacs_utils PUBLIC ${GROMACS_SOURCES}/src)
    target_link_libraries(pteros_gromacs_utils PRIVATE gromacs)
else()
    target_link_libraries(pteros_gromacs_utils PRIVATE xdrfile)
endif()
