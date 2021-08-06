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
    # These directories and defines are exposed as public!
    target_compile_definitions(pteros_gromacs_utils PUBLIC USE_GROMACS)
    target_include_directories(pteros_gromacs_utils PUBLIC
        ${GROMACS_INCLUDE_DIRECTORIS}
        ${CMAKE_BINARY_DIR}/src/core/gromacs_utils/  # For generated gromacs_version_info.h
    )
    target_link_libraries(pteros_gromacs_utils PUBLIC ${GROMACS_LIBRARIES})

    # In case of inner Gromacs build we need to set the dependency
    if(TARGET Gromacs_external)
        add_dependencies(pteros_gromacs_utils Gromacs_external)
    endif()

else()
    target_link_libraries(pteros_gromacs_utils PRIVATE xdrfile)
endif()
