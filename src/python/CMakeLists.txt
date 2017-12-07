if(WITH_PYTHON)    
    install(DIRECTORY pteros_analysis_plugins DESTINATION ${PY_INST_DIR} USE_SOURCE_PERMISSIONS)
    install(PROGRAMS
            scripts/pteros_analysis.py
            scripts/tpr2pteros.py
            DESTINATION bin)

    # Install wrapper for main extension modules
    install(DIRECTORY pteros DESTINATION ${PY_INST_DIR})

    add_subdirectory(bindings)

endif()
