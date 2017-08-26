if(WITH_PYTHON)    
    install(DIRECTORY pteros_analysis_plugins DESTINATION ${PY_INST_DIR} USE_SOURCE_PERMISSIONS)
    install(PROGRAMS
            scripts/pteros_analysis.py
            scripts/tpr2pteros.py
            DESTINATION bin)
    # Install wrapper for main extension module
    install(PROGRAMS scripts/pteros.py DESTINATION ${PY_INST_DIR})


    add_subdirectory(pybind11_src)
    INCLUDE_DIRECTORIES(pybind11_src/include)
    add_subdirectory(bindings11)

endif()
