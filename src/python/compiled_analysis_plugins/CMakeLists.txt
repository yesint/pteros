SET(PLUGINS    
    box
    rms
    energy
    energy_par
    #covar_matr
    #distance_matr
    secondary
    #rdf
    #example_plugin
    center
    contacts
    voronoi_par
)

IF(MAKE_STANDALONE_PLUGINS)
    message(STATUS "Building standalone plugins")

    foreach(plugin ${PLUGINS})
        add_executable(pteros_${plugin}
            ${plugin}.cpp
            ${PROJECT_SOURCE_DIR}/include/pteros/core/pteros_error.h
        )

        target_compile_definitions(pteros_${plugin} PRIVATE "STANDALONE_PLUGINS")
        target_link_libraries(pteros_${plugin} pteros pteros_analysis)

        install(TARGETS pteros_${plugin} RUNTIME DESTINATION bin/analysis)
    endforeach()
ENDIF()

IF(WITH_PYTHON)
    message(STATUS "Building plugins as python extensions")

    foreach(plugin ${PLUGINS})
        pybind11_add_module(${plugin}
            ${plugin}.cpp
            ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin.h
        )

        target_link_libraries(${plugin} PRIVATE pteros pteros_analysis ${PYTHON_LIBRARIES})

        set_target_properties(${plugin} PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/pteros_analysis_plugins"
        )
        #Installation
        install(TARGETS ${plugin}
            LIBRARY DESTINATION ${PY_INST_DIR}/pteros_analysis_plugins
            ARCHIVE DESTINATION ${PY_INST_DIR}/pteros_analysis_plugins
        )
    endforeach()
ENDIF()
