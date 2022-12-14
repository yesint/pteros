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
add_subdirectory(io)

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

#Add POWERSASA oncludes
# Add includes
target_include_directories(pteros PRIVATE ${PROJECT_SOURCE_DIR}/thirdparty/powersasa)

target_link_libraries(pteros
    PRIVATE
        pteros_interface
        pteros_io
        dssp
        sasa
        molfile_plugins
    PUBLIC
        fmt::fmt
        Eigen3::Eigen
        spdlog::spdlog
)

if(WITH_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(pteros PRIVATE OpenMP::OpenMP_CXX)
endif()

# Set target properties for exporting
target_include_directories(pteros PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

#--------------
# Installation
#--------------

install(TARGETS pteros EXPORT pterosTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
