if(WRAP_PYTHON)

    IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        SET(LIB_EXT ".pyd") # Force .pyd extension on Windows to ensure correct loading
    ELSE()
        SET(LIB_EXT ".so") # Use default on Linux
    ENDIF()
    
    # Compiled plugins:
    
    add_library(plugin_center SHARED 
    center.cpp
    )
    target_link_libraries(plugin_center pteros pteros_analysis pteros_sim ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})


    set_target_properties(plugin_center
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/pteros_analysis_plugins"
    PREFIX ""
    SUFFIX "${LIB_EXT}"
    OUTPUT_NAME center
    )

    #Installation
    install(TARGETS plugin_center
        LIBRARY DESTINATION python/pteros_analysis_plugins
        ARCHIVE DESTINATION python/pteros_analysis_plugins
    )

endif()
