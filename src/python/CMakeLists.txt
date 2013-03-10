if(WRAP_PYTHON)

    add_library(pteros_py SHARED 
    bindings.cpp

    bindings_util.h
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

    trajectory_processor_wrapper.cpp
    trajectory_processor_wrapper.h

    consumer_wrapper.cpp
    consumer_wrapper.h
    )

    target_link_libraries(pteros_py pteros pteros_analysis pteros_sim ${Boost_LIBRARIES})

    set_target_properties(pteros_py PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    PREFIX ""
    OUTPUT_NAME pteros_py)

else()
    message("Can't build python bindings!")
endif()

