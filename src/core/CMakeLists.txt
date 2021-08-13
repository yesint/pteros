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

add_subdirectory(babel_utils)
add_subdirectory(gromacs_utils)
add_subdirectory(io)

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

    ${PROJECT_SOURCE_DIR}/include/pteros/core/atom_handler.h
    atom_handler.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/system.h
    system.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/selection.h
    selection.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/grid.h
    grid.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/periodic_box.h
    periodic_box.cpp

    #DSSP wrapper
    pteros_dssp_wrapper.cpp
    pteros_dssp_wrapper.h
)

# Include sources in subdirectories into pteros target
add_subdirectory(distance_search)
add_subdirectory(selection_parser)

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
    message(NOTICE "POWERSASA code is used! Licence restrictions are described here: src/thirdparty/sasa/LICENSE")
endif()

target_link_libraries(pteros
    PRIVATE
        pteros_io
        dssp
        sasa
        molfile_plugins
    PUBLIC
        Eigen3::Eigen
        spdlog::spdlog
        #fmt::fmt
)

if(WITH_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(pteros PRIVATE OpenMP::OpenMP_CXX)
endif()

# Set target properties for exporting
target_include_directories(pteros PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Precompile some headers
target_precompile_headers(pteros PRIVATE
    <string>
    <memory>
    <map>
    <vector>
    <functional>
    <Eigen/Core>
    ${PROJECT_SOURCE_DIR}/include/pteros/core/selection.h
    ${PROJECT_SOURCE_DIR}/include/pteros/core/system.h
)

#--------------
# Installation
#--------------

install(TARGETS pteros EXPORT pterosTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
