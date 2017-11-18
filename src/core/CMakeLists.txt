#------------------------
# SASA code
#------------------------

if(WITH_POWERSASA)
    # Set definition for conditional compilation
    add_definitions(-DUSE_POWERSASA)
    # Add includes
    include_directories(${PROJECT_SOURCE_DIR}/thirdparty/sasa)
    # Make list of SASA files
    SET(SASA_FILES
        ${PROJECT_SOURCE_DIR}/thirdparty/sasa/power_diagram.h
        ${PROJECT_SOURCE_DIR}/thirdparty/sasa/power_sasa.h
    )

    message("------------------------------")
    message("POWERSASA code is used! \n This implies severe license restrictions! \n Read the comment in src/thirdparty/sasa/LICENSE carefully. \n If you don't agree with this license please use flag -DPOWERSASA=\"OFF\" !")
    message("------------------------------")

else()
    SET(SASA_FILES "")
endif()


#---------------------
# Main pteros library
#---------------------

add_library(pteros SHARED
    # General
    ${PROJECT_SOURCE_DIR}/include/pteros/core/logging.h
    logging.cpp
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

    ${PROJECT_SOURCE_DIR}/include/pteros/core/grid.h
    grid.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/periodic_box.h
    periodic_box.cpp

    # Distance search
    ${PROJECT_SOURCE_DIR}/include/pteros/core/distance_search.h
    distance_search/distance_search.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/distance_search_within.h
    distance_search/distance_search_within.cpp

    distance_search/atomic_wrapper.h
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

    #SASA (will be empty if not used)
    ${SASA_FILES}

    #DSSP wrapper
    ${PROJECT_SOURCE_DIR}/src/core/pteros_dssp_wrapper.cpp
    ${PROJECT_SOURCE_DIR}/src/core/pteros_dssp_wrapper.h
)

target_link_libraries(pteros tng_io dssp molfile_plugins xdrfile)

#--------------
# Installation
#--------------

install(TARGETS pteros    
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
