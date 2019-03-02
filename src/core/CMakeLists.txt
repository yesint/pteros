# Obtain version from git
find_package(Git)
find_file(GITDIR NAMES .git PATHS ${PROJECT_SOURCE_DIR} NO_DEFAULT_PATH)
if(GIT_FOUND AND GITDIR)
    SET(_git_head_dependence "${GITDIR}/logs/HEAD")
else()
    SET(_git_head_dependence "")
endif()
add_custom_command(OUTPUT "${PROJECT_SOURCE_DIR}/src/core/version.cpp"
    COMMAND ${CMAKE_COMMAND}
        -Dwork_dir=${PROJECT_SOURCE_DIR}
        -P "${PROJECT_SOURCE_DIR}/cmake/gitversion.cmake"
    DEPENDS ${_git_head_dependence}
    VERBATIM
)


#------------------------
# SASA code
#------------------------

add_subdirectory(io)
add_subdirectory(selection_parser)
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
    ${PROJECT_SOURCE_DIR}/include/pteros/core/version.h
    version.cpp # generated with git revision if git present

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
    selection.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/grid.h
    grid.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/periodic_box.h
    periodic_box.cpp

    #SASA (will be empty if not used)
    ${SASA_FILES}

    #DSSP wrapper
    pteros_dssp_wrapper.cpp
    pteros_dssp_wrapper.h
)

# Add path to molfile plugins directory to use VMD periodic table routines
target_include_directories(pteros PRIVATE ${PROJECT_SOURCE_DIR}/thirdparty/molfile_plugins)

# For selection parser includes
target_include_directories(pteros PRIVATE selection_parser)

target_link_libraries(pteros pteros_io pteros_selection_parser pteros_distance_search dssp sasa)

#--------------
# Installation
#--------------

install(TARGETS pteros    
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
