if(WITH_PYTHON)
    file(GLOB scripts_to_install "scripts/*.py")
    install(PROGRAMS
            ${scripts_to_install}
            DESTINATION bin/pteros)

    # Install main python module directory
    install(DIRECTORY pteros DESTINATION ${PY_INST_DIR})

    add_subdirectory(bindings)
endif()
