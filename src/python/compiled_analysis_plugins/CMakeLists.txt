
SET(PLUGINS
    center
    box
    rms
    interaction_energy
    covar_matr
    distance_matr
    secondary
    rdf
)

IF(STANDALONE_PLUGINS)

    message("-----------------------------------------")
    message("      Building standalone plugins        ")
    message("-----------------------------------------")


    # Conditional compilation definition
    #add_definitions(-DSTANDALONE_PLUGINS)

    foreach(plugin ${PLUGINS})
        add_executable(pteros_${plugin}
            ${plugin}.cpp
            ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin_base.h
            ${PROJECT_SOURCE_DIR}/include/pteros/core/pteros_error.h
        )
        #set_source_files_properties(${plugin}.cpp  PROPERTIES COMPILE_DEFINITIONS "STANDALONE_PLUGINS")
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
        add_library(plugin_${plugin} SHARED
            ${plugin}.cpp
            ${PROJECT_SOURCE_DIR}/include/pteros/python/bindings_util.h
            ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin_base.h
            ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin.h
        )
        target_link_libraries(plugin_${plugin} pteros pteros_analysis ${PYTHON_LIBRARIES} ${BOOST_PYTHON_LIB})
        set_target_properties(plugin_${plugin}
            PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/pteros_analysis_plugins"
            PREFIX "" SUFFIX "${SHARED_LIB_EXT}"
            OUTPUT_NAME ${plugin}
        )
        #Installation
        install(TARGETS plugin_${plugin}
            LIBRARY DESTINATION python/pteros_analysis_plugins
            ARCHIVE DESTINATION python/pteros_analysis_plugins
        )
    endforeach()

ENDIF()
