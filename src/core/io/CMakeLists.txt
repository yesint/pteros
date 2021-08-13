#---------------------
# pteros IO library
#---------------------

add_library(pteros_io STATIC
    ${PROJECT_SOURCE_DIR}/include/pteros/core/file_handler.h
    file_handler.cpp
    system_builder.h
    system_builder.cpp
    vmd_molfile_plugin_wrapper.h
    vmd_molfile_plugin_wrapper.cpp
    pdb_file.h
    pdb_file.cpp
    dcd_file.h
    dcd_file.cpp
    gro_file.h
    gro_file.cpp    
    xyz_file.h
    xyz_file.cpp
    trr_file.h
    trr_file.cpp
    xtc_file.h
    xtc_file.cpp
)

if(WITH_TNG)
    target_sources(pteros_io PRIVATE
        tng_file.h
        tng_file.cpp)
    target_compile_definitions(pteros_io PRIVATE USE_TNGIO)
    target_link_libraries(pteros_io PRIVATE ${TNG_LIBRARIES})
    target_include_directories(pteros_io PRIVATE ${TNG_INCLUDE_DIR})

    if(TARGET TNG_external)
        add_dependencies(pteros_io TNG_external)
    endif()
endif()

if(WITH_OPENBABEL)
    target_sources(pteros_io PRIVATE
        babel_wrapper.h
        babel_wrapper.cpp
        pdbqt_file.h
        pdbqt_file.cpp
        mol2_file.h
        mol2_file.cpp)
    target_compile_definitions(pteros_io PRIVATE USE_OPENBABEL)
    target_include_directories(pteros_io PRIVATE ${OPENBABEL2_INCLUDE_DIR} ${OPENBABEL3_INCLUDE_DIR})
    target_link_libraries(pteros_io PRIVATE pteros_babel_utils)
endif()

if(WITH_GROMACS)
    # All definitions and includes are given by pteros_gromacs_utils
    target_sources(pteros_io PRIVATE
        tpr_file.h
        tpr_file.cpp)
endif()

# Unconditional libs
target_link_libraries(pteros_io PRIVATE
    pteros_gromacs_utils  # Brings all Gromacs defines and directories
    molfile_plugins
    xdrfile
    spdlog::spdlog
    Eigen3::Eigen
    #fmt::fmt
)
