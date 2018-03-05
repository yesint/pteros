#------------------------
# SASA code
#------------------------

add_subdirectory(io)
add_subdirectory(distance_search)

if(WITH_POWERSASA)
    # Set definition for conditional compilation
    add_definitions(-DUSE_POWERSASA)
    # Add includes
    include_directories(${PROJECT_SOURCE_DIR}/thirdparty/powersasa)
    # Make list of SASA files
    SET(SASA_FILES
        ${PROJECT_SOURCE_DIR}/thirdparty/powersasa/power_diagram.h
        ${PROJECT_SOURCE_DIR}/thirdparty/powersasa/power_sasa.h
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
    ${PROJECT_SOURCE_DIR}/include/pteros/core/atom_proxy.h
    atom_proxy.cpp
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

    #SASA (will be empty if not used)
    ${SASA_FILES}

    #DSSP wrapper
    ${PROJECT_SOURCE_DIR}/src/core/pteros_dssp_wrapper.cpp
    ${PROJECT_SOURCE_DIR}/src/core/pteros_dssp_wrapper.h
)

# Add path to molfile plugins directory to use VMD periodic table routines
target_include_directories(pteros PRIVATE ${PROJECT_SOURCE_DIR}/thirdparty/molfile_plugins)

target_link_libraries(pteros pteros_io pteros_distance_search dssp sasa)

#--------------
# Installation
#--------------

install(TARGETS pteros    
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
