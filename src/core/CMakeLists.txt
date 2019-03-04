# Obtain version from git and generate version.cpp
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
    DEPENDS version.cpp.in  ${_git_head_dependence}
    VERBATIM
)

add_subdirectory(io)
add_subdirectory(selection_parser)
add_subdirectory(distance_search)
add_subdirectory(babel_utils)

#-------------------------------------------------------------------------------
# service interface library to provide includes for other targets at build time
#-------------------------------------------------------------------------------

add_library(pteros_interface INTERFACE)
target_include_directories(pteros_interface INTERFACE ${PROJECT_SOURCE_DIR}/include/)

#---------------------
# pteros core library
#---------------------

add_library(pteros SHARED
    # General    
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

#Add SASA code

if(WITH_POWERSASA)
    # Set definition for conditional compilation
    target_compile_definitions(pteros PRIVATE USE_POWERSASA)
    # Add includes
    target_include_directories(pteros PRIVATE ${PROJECT_SOURCE_DIR}/thirdparty/powersasa)
    # sources
    target_sources(pteros PRIVATE
        ${PROJECT_SOURCE_DIR}/thirdparty/powersasa/power_diagram.h
        ${PROJECT_SOURCE_DIR}/thirdparty/powersasa/power_sasa.h
    )
    message(STATUS "POWERSASA code is used! Licence restrictions are described here: src/thirdparty/sasa/LICENSE")
endif()

target_link_libraries(pteros
    PRIVATE
        pteros_io
        pteros_selection_parser
        pteros_distance_search
        dssp
        sasa
        molfile_plugins
        Boost::boost Boost::system
    PUBLIC
        Eigen3::Eigen
        spdlog::spdlog
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(pteros PRIVATE OpenMP::OpenMP_CXX)
endif()

# Set target properties for exporting

target_include_directories(pteros PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

set_property(TARGET pteros
    PROPERTY CXX_STANDARD 14
    PROPERTY CXX_STANDARD_REQUIRED ON
)

#--------------
# Installation
#--------------

install(TARGETS pteros EXPORT pterosTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
