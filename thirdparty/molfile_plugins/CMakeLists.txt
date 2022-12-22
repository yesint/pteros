#---------------------
# VMD molfile plugins
#---------------------

#---------------------------
# Build machinery
#---------------------------

SET(MOLFILE_PLUGINS pdbplugin dcdplugin xyzplugin)
if(WITH_TNG)
    list(APPEND MOLFILE_PLUGINS tngplugin)
    set_source_files_properties(register_plugins.c PROPERTIES COMPILE_DEFINITIONS USE_TNG)
endif()

SET(MOLFILE_PLUGINS_FILE_LIST "")

foreach(PLUGIN ${MOLFILE_PLUGINS})
    file(GLOB F ${PLUGIN}.*)
    set_source_files_properties(${F} PROPERTIES COMPILE_DEFINITIONS VMDPLUGIN=${PLUGIN})
    list(APPEND MOLFILE_PLUGINS_FILE_LIST ${F})
endforeach(PLUGIN)

add_library(molfile_plugins STATIC
    register_plugins.c
    ${MOLFILE_PLUGINS_FILE_LIST}
)

target_include_directories(molfile_plugins PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(WITH_TNG)
    target_link_directories(molfile_plugins PUBLIC tng_interface)
    # For some reason includes from tng_interface are not used...
    target_include_directories(molfile_plugins PRIVATE ${TNG_INCLUDE_DIR})
endif()
