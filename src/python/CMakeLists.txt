if(WITH_PYTHON)    
    install(DIRECTORY pteros_analysis_plugins DESTINATION ${PY_INST_DIR} USE_SOURCE_PERMISSIONS)
    install(PROGRAMS
            scripts/pteros_analysis.py
            scripts/tpr2pteros.py
            DESTINATION bin)
    # Install wrapper for main extension module
    install(PROGRAMS scripts/pteros.py DESTINATION ${PY_INST_DIR})


    add_subdirectory(pybind11)
    INCLUDE_DIRECTORIES(pybind11/include)
    pybind11_add_module(pteros11 bindings.cpp)

    target_link_libraries(pteros11 PRIVATE pteros pteros_analysis)

    set_target_properties(pteros11 PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    )

    #Installation
    install(TARGETS pteros11
        LIBRARY DESTINATION ${PY_INST_DIR}
        ARCHIVE DESTINATION ${PY_INST_DIR}
    )


endif()
