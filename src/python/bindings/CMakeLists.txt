IF(WITH_PYTHON)

    ###########################
    # Main pteros python module
    ###########################

    add_library(pteros_py SHARED 
        bindings.cpp

        ${PROJECT_SOURCE_DIR}/include/pteros/python/bindings_util.h
        ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin_base.h
        ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin.h

        bindings_system.cpp
        bindings_system.h
        bindings_selection.cpp
        bindings_selection.h
        bindings_utilities.cpp
        bindings_utilities.h
        bindings_frame_info.cpp
        bindings_frame_info.h
        bindings_trajectory_processor.cpp
        bindings_trajectory_processor.h
        bindings_atom.cpp
        bindings_atom.h
        bindings_atom_proxy.cpp
        bindings_atom_proxy.h
        bindings_distance_search.cpp
        bindings_distance_search.h
        bindings_periodic_box.cpp
        bindings_periodic_box.h
        bindings_frame.cpp
        bindings_frame.h
        bindings_energy_components.cpp
        bindings_energy_components.h

        trajectory_processor_wrapper.cpp
        trajectory_processor_wrapper.h

        consumer_wrapper.cpp
        consumer_wrapper.h

        bindings_options.cpp
        bindings_options.h
    )

    target_link_libraries(pteros_py pteros pteros_analysis ${PYTHON_LIBRARIES} ${BOOST_PYTHON_LIB})

    set_target_properties(pteros_py PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
        PREFIX ""
        SUFFIX "${SHARED_LIB_EXT}"
        OUTPUT_NAME _pteros
    )

    #Installation
    install(TARGETS pteros_py
        LIBRARY DESTINATION ${PY_INST_DIR}
        ARCHIVE DESTINATION ${PY_INST_DIR}
    )

ENDIF()

