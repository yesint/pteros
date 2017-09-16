#---------------------
# VMD molfile plugins
#---------------------

#---------------------------
# Build machinery
#---------------------------

SET(MOLFILE_PLUGINS pdbplugin tngplugin dcdplugin mol2plugin xyzplugin)
SET(MOLFILE_PLUGINS_FILE_LIST "")

foreach(PLUGIN ${MOLFILE_PLUGINS})
    file(GLOB F ${PLUGIN}.*)
    set_source_files_properties(${F} PROPERTIES COMPILE_DEFINITIONS VMDPLUGIN=${PLUGIN})
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

