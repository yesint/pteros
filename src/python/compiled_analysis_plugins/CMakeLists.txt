if(WRAP_PYTHON)

    IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        SET(LIB_EXT ".pyd") # Force .pyd extension on Windows to ensure correct loading
    ELSE()
        SET(LIB_EXT ".so") # Use default on Linux
    ENDIF()
    
    link_libraries(pteros pteros_analysis ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

    # Compiled plugins:
    IF(NO_CPP11)
        SET(PLUGINS center box rmsd interaction_energy covar_matr distance_matr)
    ELSE()
        SET(PLUGINS center box rmsd interaction_energy covar_matr distance_matr secondary)
    ENDIF()

    foreach(plugin ${PLUGINS})
        add_library(plugin_${plugin} SHARED ${plugin}.cpp)
        set_target_properties(plugin_${plugin}
        PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/pteros_analysis_plugins"
        PREFIX "" SUFFIX "${LIB_EXT}"
        OUTPUT_NAME ${plugin} )
        #Installation
        install(TARGETS plugin_${plugin}
            LIBRARY DESTINATION python/pteros_analysis_plugins
            ARCHIVE DESTINATION python/pteros_analysis_plugins
        )
    endforeach()


endif()
