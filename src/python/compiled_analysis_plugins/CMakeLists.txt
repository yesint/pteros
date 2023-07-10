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
    charge_distrib
)

IF(MAKE_STANDALONE_PLUGINS)
    message(STATUS "Building standalone plugins")

    foreach(plugin ${PLUGINS})
        add_executable(
            pteros_${plugin}
            ${plugin}.cpp
        )

        target_compile_definitions(pteros_${plugin} PRIVATE "STANDALONE_PLUGINS")
        target_link_libraries(pteros_${plugin} PRIVATE pteros pteros_analysis)

        install(TARGETS pteros_${plugin} RUNTIME DESTINATION bin/analysis)
    endforeach()

    # Extra link dependency for voronoi plugin
    target_link_libraries(pteros_voronoi_par PRIVATE pteros_voronoi_packing)
ENDIF()

IF(WITH_PYTHON)
    message(STATUS "Building plugins as python extensions")

    foreach(plugin ${PLUGINS})
        pybind11_add_module(
            ${plugin}
            ${plugin}.cpp
        )

        target_link_libraries(${plugin} PRIVATE pteros pteros_analysis ${PYTHON_LIBRARIES})

        #Installation
        install(TARGETS ${plugin}
            LIBRARY DESTINATION ${PY_INST_DIR}/pteros/plugins
            ARCHIVE DESTINATION ${PY_INST_DIR}/pteros/plugins
        )
    endforeach()

    # Extra link dependency for voronoi plugin
    target_link_libraries(voronoi_par PRIVATE pteros_voronoi_packing)
ENDIF()
