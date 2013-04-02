if(WRAP_PYTHON)
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        SET(LIB_EXT ".pyd") # Force .pyd extension on Windows to ensure correct loading
    ELSE()
        SET(LIB_EXT ".so") # Use default on Linux
    ENDIF()

    add_library(pteros_py SHARED 
    bindings.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/python/bindings_util.h
    ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin_base.h
    ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin.h
    bindings_system.cpp
    bindings_system.h
    bindings_selection.cpp
    bindings_selection.h
    bindings_options_tree.cpp
    bindings_options_tree.h
    bindings_frame_info.cpp
    bindings_frame_info.h
    bindings_trajectory_processor.cpp
    bindings_trajectory_processor.h
    bindings_atom.cpp
    bindings_atom.h

    trajectory_processor_wrapper.cpp
    trajectory_processor_wrapper.h

    consumer_wrapper.cpp
    consumer_wrapper.h
    )

    target_link_libraries(pteros_py pteros pteros_analysis pteros_sim ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

    set_target_properties(pteros_py PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    PREFIX ""
    SUFFIX "${LIB_EXT}"
    OUTPUT_NAME pteros)

    #Installation
    install(TARGETS pteros_py
        LIBRARY DESTINATION python
        ARCHIVE DESTINATION python
    )


else()
    message("Python bindings disabled!")
endif()
