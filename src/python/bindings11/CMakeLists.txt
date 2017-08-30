IF(WITH_PYTHON)

    #############################
    # Main pteros11 python module
    #############################
    SET(BIND_FILES
        bindings.cpp
        bindings_util.h
        bindings_atom.cpp
        bindings_periodic_box.cpp
        bindings_system.cpp
        bindings_selection.cpp
        bindings_frame.cpp
        bindings_distance_search.cpp
        bindings_options.cpp
        bindings_trajectory_reader.cpp
    )


    pybind11_add_module(pteros11 ${BIND_FILES})

    target_link_libraries(pteros11 PRIVATE pteros pteros_analysis)

    set_target_properties(pteros11 PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    )

    #Installation
    install(TARGETS pteros11
        LIBRARY DESTINATION ${PY_INST_DIR}
        ARCHIVE DESTINATION ${PY_INST_DIR}
    )

ENDIF()

