if(WITH_OPENBABEL AND ${OPENBABEL2_FOUND})

    add_library(pteros_babel_utils STATIC
        babel_utils.cpp
        babel_utils.h
    )

    target_link_libraries(pteros_babel_utils
        PUBLIC
            ${OPENBABEL2_LIBRARIES}
            Eigen3::Eigen
            pteros_core
    )

    target_include_directories(pteros_babel_utils
        PUBLIC
            ${OPENBABEL2_INCLUDE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}
    )

endif()
