
SET(PLUGINS    
    box
    rms
    energy
    #covar_matr
    #distance_matr
    #secondary
    #rdf
    #example_plugin
    center
    contacts
)

IF(STANDALONE_PLUGINS)

    message("-----------------------------------------")
    message("      Building standalone plugins        ")
    message("-----------------------------------------")

    foreach(plugin ${PLUGINS})
        add_executable(pteros_${plugin}
            ${plugin}.cpp            
            ${PROJECT_SOURCE_DIR}/include/pteros/core/pteros_error.h
        )

        # Set compile definitions for this target
        set_target_properties(pteros_${plugin} PROPERTIES COMPILE_DEFINITIONS "STANDALONE_PLUGINS")

        target_link_libraries(pteros_${plugin} pteros pteros_analysis)

        install(TARGETS pteros_${plugin} RUNTIME DESTINATION bin/analysis)
    endforeach()

ENDIF()

IF(WITH_PYTHON)

    message("-----------------------------------------")
    message("  Building plugins as python extensions  ")
    message("-----------------------------------------")

    foreach(plugin ${PLUGINS})
        pybind11_add_module(${plugin}
            ${plugin}.cpp
            ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin.h
        )

        target_link_libraries(${plugin} PRIVATE pteros pteros_analysis ${PYTHON_LIBRARIES})

        set_target_properties(${plugin}
            PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/pteros_analysis_plugins"
        )
        #Installation
        install(TARGETS ${plugin}
            LIBRARY DESTINATION ${PY_INST_DIR}/pteros_analysis_plugins
            ARCHIVE DESTINATION ${PY_INST_DIR}/pteros_analysis_plugins
        )
    endforeach()

ENDIF()
