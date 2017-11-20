IF(WITH_PYTHON)

    #############################
    # Main pteros python module
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
        bindings_membrane.cpp
        ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin.h
    )


    pybind11_add_module(_pteros ${BIND_FILES})

    target_link_libraries(_pteros PRIVATE pteros pteros_analysis pteros_membrane)

    set_target_properties(_pteros PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    )

    #Installation
    install(TARGETS _pteros
        LIBRARY DESTINATION ${PY_INST_DIR}
        ARCHIVE DESTINATION ${PY_INST_DIR}
    )

ENDIF()

