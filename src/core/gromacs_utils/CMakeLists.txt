add_library(pteros_gromacs_utils STATIC
    gromacs_utils.cpp
    gromacs_utils.h
)

target_link_libraries(pteros_gromacs_utils
    PUBLIC
        Eigen3::Eigen
        pteros_interface
)

target_include_directories(pteros_gromacs_utils
    PUBLIC
        ${GROMACS_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

if(WITH_GROMACS AND gromacs_FOUND)
    target_compile_definitions(pteros_gromacs_utils PRIVATE USE_GROMACS)
    target_link_libraries(pteros_gromacs_utils PRIVATE ${GROMACS_LIBRARIES})
else()
    target_link_libraries(pteros_gromacs_utils PRIVATE xdrfile)
endif()
