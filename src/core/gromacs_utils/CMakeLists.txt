if(WITH_GROMACS AND gromacs_FOUND)

    add_library(pteros_gromacs_utils STATIC
        gromacs_utils.cpp
        gromacs_utils.h
    )

    target_link_libraries(pteros_gromacs_utils
        PUBLIC
            ${GROMACS_LIBRARIES}
            Eigen3::Eigen
            pteros_interface
    )

    target_include_directories(pteros_gromacs_utils
        PUBLIC
            ${GROMACS_INCLUDE_DIRS}
            ${CMAKE_CURRENT_SOURCE_DIR}
    )

endif()
