#---------------------
# pteros IO library
#---------------------

add_library(pteros_io STATIC
    ${PROJECT_SOURCE_DIR}/include/pteros/core/mol_file.h
    mol_file.cpp
    vmd_molfile_plugin_wrapper.h
    vmd_molfile_plugin_wrapper.cpp

    pdb_file.h
    pdb_file.cpp
    dcd_file.h
    dcd_file.cpp
    gro_file.h
    gro_file.cpp
    gromacs_trajectory_file.h
    gromacs_trajectory_file.cpp
    trr_file.h
    trr_file.cpp
    xtc_file.h
    xtc_file.cpp    
    tng_file.h
    tng_file.cpp
    mol2_file.h
    mol2_file.cpp
    xyz_file.h
    xyz_file.cpp
    tpr_file.h
    tpr_file.cpp    
)

if(WITH_OPENBABEL AND OPENBABEL2_FOUND)
    target_sources(pteros_io PRIVATE
        babel_wrapper.h
        babel_wrapper.cpp
        pdbqt_file.h
        pdbqt_file.cpp)

    target_compile_definitions(pteros_io PRIVATE USE_OPENBABEL)
    target_include_directories(pteros_io PRIVATE ${OPENBABEL2_INCLUDE_DIR})
    target_link_libraries(pteros_io
        PRIVATE
            ${OPENBABEL2_LIBRARIES}
            pteros_babel_utils
    )
endif()

target_link_libraries(pteros_io PRIVATE
    molfile_plugins
    xdrfile
    tng_io
    spdlog::spdlog
    Boost::boost
    Eigen3::Eigen)

if(${gromacs_FOUND})
    # Public to propagate to pteros core
    target_compile_definitions(pteros_io PRIVATE USE_GROMACS)
    target_include_directories(pteros_io PRIVATE ${GROMACS_INCLUDE_DIRS})
    target_link_libraries(pteros_io PRIVATE ${GROMACS_LIBRARIES})
    # If we have version < 2018 set flag for old API
    if(${gromacs_VERSION_MAJOR} LESS 2018)
        target_compile_definitions(pteros_io PRIVATE GROMACS_OLD)
        message(STATUS "Gromacs version ${gromacs_VERSION_MAJOR} thus using old gmx_mtop_t_to_t_topology API")
    endif()
endif()
