# Standalone plugins are incompatible with python plugins, so if they are enabled
# no python plugins are built
IF(STANDALONE_PLUGINS)

    message("-----------------------------------------")
    message("(!WARNING!) Building standalone plugins!")
    message("This disables building of python compiled plugins!")

    link_libraries(pteros pteros_analysis ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

    SET(CMAKE_CXX_FLAGS "${OLD_CXX_FLAGS} -DSTANDALONE_PLUGINS")

    foreach(plugin ${PLUGINS})
        add_executable(pteros_${plugin} ${plugin}.cpp)
        install(TARGETS pteros_${plugin}
            RUNTIME DESTINATION bin/analysis
        )
    endforeach()

ELSE()

    if(WRAP_PYTHON)
        IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
            SET(LIB_EXT ".pyd") # Force .pyd extension on Windows to ensure correct loading
        ELSE()
            SET(LIB_EXT ".so") # Use default on Linux
        ENDIF()

        # List of compiled plugins:
        IF(NO_CPP11)
            SET(PLUGINS center box rms interaction_energy covar_matr distance_matr)
        ELSE()
            SET(PLUGINS center box rms interaction_energy covar_matr distance_matr secondary)
        ENDIF()

        foreach(plugin ${PLUGINS})
            add_library(plugin_${plugin} SHARED
                ${plugin}.cpp
                ${PROJECT_SOURCE_DIR}/include/pteros/python/bindings_util.h
                ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin_base.h
                ${PROJECT_SOURCE_DIR}/include/pteros/python/compiled_plugin.h
            )
            target_link_libraries(plugin_${plugin} pteros pteros_analysis ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
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

ENDIF()
