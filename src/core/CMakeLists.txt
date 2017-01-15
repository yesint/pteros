# Build separate tng_io library (will be linked later)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/core/tng_io)

#------------------------
# SASA code
#------------------------

if(POWERSASA)
    # Set definition for conditional compilation
    add_definitions(-DUSE_POWERSASA)
    # Make list of SASA files
    SET(SASA_FILES
        ${PROJECT_SOURCE_DIR}/src/core/sasa/power_diagram.h
        ${PROJECT_SOURCE_DIR}/src/core/sasa/power_sasa.h
    )

    message("------------------------------")
    message("POWERSASA code is used! \n This implies severe license restrictions! \n Read the comment in src/core/sasa/power_sasa.h carefully. \n If you don't agree with this license please use flag -DPOWERSASA=\"OFF\" !")

else()
    SET(SASA_FILES "")
endif()

#---------------------------------------------
# Set locations of headers for 3rd party libs
#---------------------------------------------

# set includes for SASA
include_directories(${PROJECT_SOURCE_DIR}/src/core/sasa)
# Headers of xdrfile are stored together with sources
include_directories(${PROJECT_SOURCE_DIR}/src/core/xdrfile)
# The same is for molfile plugins
include_directories(${PROJECT_SOURCE_DIR}/src/core/molfile_plugins)
# The same is for DSSP
include_directories(${PROJECT_SOURCE_DIR}/src/core/dssp)
# For tng
include_directories(${PROJECT_SOURCE_DIR}/src/core/tng_io/include/tng)

#---------------------------
# Molfile plugins machinery
#---------------------------

SET(MOLFILE_PLUGINS_DIR ${PROJECT_SOURCE_DIR}/src/core/molfile_plugins)

SET(MOLFILE_PLUGINS pdbplugin tngplugin dcdplugin mol2plugin xyzplugin)
SET(MOLFILE_PLUGINS_FILE_LIST "")

foreach(PLUGIN ${MOLFILE_PLUGINS})
    file(GLOB F ${MOLFILE_PLUGINS_DIR}/${PLUGIN}.*)    
    set_source_files_properties(${F} PROPERTIES COMPILE_DEFINITIONS VMDPLUGIN=${PLUGIN})
    list(APPEND MOLFILE_PLUGINS_FILE_LIST ${F})
endforeach(PLUGIN)


#---------------------
# Main pteros library
#---------------------

add_library(pteros SHARED
    # General
    ${PROJECT_SOURCE_DIR}/include/pteros/core/utilities.h
    utilities.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/pteros.h
    ${PROJECT_SOURCE_DIR}/include/pteros/core/typedefs.h
    ${PROJECT_SOURCE_DIR}/include/pteros/core/pteros_error.h
    ${PROJECT_SOURCE_DIR}/include/pteros/core/atom.h
    ${PROJECT_SOURCE_DIR}/include/pteros/core/force_field.h
    force_field.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/core/system.h
    system.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/core/selection.h
    selection_macro.h    
    selection_parser.h
    selection_grammar.h
    selection_parser.cpp
    selection.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/core/gnm.h
    gnm.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/core/grid.h
    grid.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/periodic_box.h
    periodic_box.cpp

    # Distance search
    ${PROJECT_SOURCE_DIR}/include/pteros/core/distance_search.h
    distance_search/distance_search.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/distance_search_within.h
    distance_search/distance_search_within.cpp

    distance_search/search_utils.h
    distance_search/search_utils.cpp

    distance_search/distance_search_base.h
    distance_search/distance_search_base.cpp

    distance_search/distance_search_within_base.h
    distance_search/distance_search_within_base.cpp

    distance_search/distance_search_within_sel.h
    distance_search/distance_search_within_sel.cpp

    distance_search/distance_search_contacts.h
    distance_search/distance_search_contacts.cpp

    distance_search/distance_search_contacts_1sel.h
    distance_search/distance_search_contacts_1sel.cpp

    distance_search/distance_search_contacts_2sel.h
    distance_search/distance_search_contacts_2sel.cpp


    # IO
    ${PROJECT_SOURCE_DIR}/include/pteros/core/mol_file.h
    io/mol_file.cpp
    io/vmd_molfile_plugin_wrapper.h
    io/vmd_molfile_plugin_wrapper.cpp
    io/pdb_file.h
    io/pdb_file.cpp
    io/dcd_file.h
    io/dcd_file.cpp
    io/gro_file.h
    io/gro_file.cpp
    io/gromacs_trajectory_file.h
    io/gromacs_trajectory_file.cpp
    io/trr_file.h
    io/trr_file.cpp
    io/xtc_file.h
    io/xtc_file.cpp
    io/pttop_file.h
    io/pttop_file.cpp
    io/tng_file.h
    io/tng_file.cpp
    io/mol2_file.h
    io/mol2_file.cpp
    io/xyz_file.h
    io/xyz_file.cpp

    # xdrfile
    ${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile.c
    ${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile_trr.c
    ${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile_xtc.c
    ${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile.h
    ${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile_trr.h
    ${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile_xtc.h

    # molfile_plugins
    # required headers
    ${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/molfile_plugin.h
    ${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/periodic_table.h
    ${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/vmdplugin.h
    ${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/largefiles.h
    ${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/endianswap.h    
    ${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/readpdb.h
    ${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/fastio.h
    # plugins themselves
    ${MOLFILE_PLUGINS_FILE_LIST}

    #SASA (will be empty if not used)
    ${SASA_FILES}

    #DSSP
    ${PROJECT_SOURCE_DIR}/src/core/dssp/pteros_dssp_wrapper.cpp
    ${PROJECT_SOURCE_DIR}/src/core/dssp/pteros_dssp_wrapper.h    
    ${PROJECT_SOURCE_DIR}/src/core/dssp/align-2d.h
    ${PROJECT_SOURCE_DIR}/src/core/dssp/matrix.h
    ${PROJECT_SOURCE_DIR}/src/core/dssp/dssp.h
    ${PROJECT_SOURCE_DIR}/src/core/dssp/mas.h
    ${PROJECT_SOURCE_DIR}/src/core/dssp/primitives-3d.h
    ${PROJECT_SOURCE_DIR}/src/core/dssp/structure.h
    ${PROJECT_SOURCE_DIR}/src/core/dssp/utils.h
    ${PROJECT_SOURCE_DIR}/src/core/dssp/dssp.cpp
    ${PROJECT_SOURCE_DIR}/src/core/dssp/mas.cpp
    ${PROJECT_SOURCE_DIR}/src/core/dssp/primitives-3d.cpp
    ${PROJECT_SOURCE_DIR}/src/core/dssp/structure.cpp
    ${PROJECT_SOURCE_DIR}/src/core/dssp/utils.cpp    
)

target_link_libraries(pteros tng_io)

#--------------
# Installation
#--------------

install(TARGETS pteros    
    RUNTIME DESTINATION lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
