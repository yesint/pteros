if(WITH_PYTHON)

    install(DIRECTORY pteros_analysis_plugins DESTINATION python USE_SOURCE_PERMISSIONS)
    install(PROGRAMS
            scripts/pteros_analysis.py
            scripts/tpr2pteros.py
            DESTINATION bin)

endif()
