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

    ####################################
    # pteros_extras python module files
    ####################################
    SET(EXTRAS_FILES
        bindings_extras.cpp
        bindings_membrane.cpp
    )


    pybind11_add_module(_pteros ${BIND_FILES})
    pybind11_add_module(_pteros_extras ${EXTRAS_FILES})

    target_link_libraries(_pteros PRIVATE pteros pteros_analysis)
    target_link_libraries(_pteros_extras PRIVATE pteros_membrane)

    #Installation
    install(TARGETS _pteros
        LIBRARY DESTINATION ${PY_INST_DIR}/pteros
        ARCHIVE DESTINATION ${PY_INST_DIR}/pteros
    )

    install(TARGETS _pteros_extras
        LIBRARY DESTINATION ${PY_INST_DIR}/pteros/extras
        ARCHIVE DESTINATION ${PY_INST_DIR}/pteros/extras
    )

ENDIF()

