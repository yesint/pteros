if(WITH_PYTHON)

    install(DIRECTORY pteros_analysis_plugins DESTINATION python USE_SOURCE_PERMISSIONS)
    install(PROGRAMS
            scripts/pteros_analysis.py
            scripts/tpr2pteros.py
            DESTINATION bin)
    # Install wrapper for main extension module
    install(PROGRAMS scripts/pteros.py DESTINATION python)
endif()
