#---------------------
# VMD molfile plugins
#---------------------

#---------------------------
# Build machinery
#---------------------------

SET(MOLFILE_PLUGINS pdbplugin dcdplugin mol2plugin xyzplugin)
if(WITH_TNGIO)
    list(APPEND MOLFILE_PLUGINS tngplugin)
endif()

SET(MOLFILE_PLUGINS_FILE_LIST "")

foreach(PLUGIN ${MOLFILE_PLUGINS})
    file(GLOB F ${PLUGIN}.*)
    set_source_files_properties(${F} PROPERTIES COMPILE_DEFINITIONS VMDPLUGIN=${PLUGIN})
    if(MINGW)
        set_source_files_properties(${F} PROPERTIES COMPILE_DEFINITIONS STATIC_PLUGIN)
    endif()
    list(APPEND MOLFILE_PLUGINS_FILE_LIST ${F})
endforeach(PLUGIN)

add_library(molfile_plugins STATIC
    # required headers
    molfile_plugin.h
    periodic_table.h
    vmdplugin.h
    largefiles.h
    endianswap.h
    readpdb.h
    fastio.h
    # plugins themselves
    ${MOLFILE_PLUGINS_FILE_LIST}
)

if(WITH_TNGIO)
    target_link_libraries(molfile_plugins PRIVATE tng_io)
endif()

target_include_directories(molfile_plugins PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
