if(WITH_OPENBABEL AND (OPENBABEL2_FOUND OR OPENBABEL3_FOUND))

    add_library(pteros_babel_utils STATIC
        babel_utils.cpp
        babel_utils.h
    )

    target_link_libraries(pteros_babel_utils
        PUBLIC
            ${OPENBABEL2_LIBRARIES} ${OPENBABEL3_LIBRARIES}
            Eigen3::Eigen
            pteros_interface
    )

    target_include_directories(pteros_babel_utils
        PUBLIC
            ${OPENBABEL2_INCLUDE_DIR} ${OPENBABEL3_INCLUDE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}
    )

    # In case of inner OpenBabel build we need to set the dependency
    if(TARGET OpenBabel_external)
        add_dependencies(pteros_babel_utils OpenBabel_external)
    endif()

endif()
