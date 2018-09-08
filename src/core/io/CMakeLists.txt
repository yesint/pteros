#---------------------
# pteros IO library
#---------------------

# Look for Gromacs
if(WITH_GROMACS)
    message("Looking for Gromacs...")
    find_package(gromacs QUIET)
    if(NOT ${gromacs_FOUND})
        # Try /usr/local
        message("   Trying /usr/local...")
        set(gromacs_DIR "/usr/local/gromacs/share/cmake/gromacs")
        find_package(gromacs QUIET)
    endif()

    if(${gromacs_FOUND})
        message("-- Found Gromacs ${gromacs_VERSION} at ${GROMACS_INCLUDE_DIRS}")
        add_definitions(-DUSE_GROMACS)
        include_directories(${GROMACS_INCLUDE_DIRS})
        link_libraries(${GROMACS_LIBRARIES})
        # If we have version < 2018 set flag for old API
        if(${gromacs_VERSION_MAJOR} LESS 2018)
            add_definitions(-DGROMACS_OLD)
            message("-- Gromacs version ${gromacs_VERSION_MAJOR} thus using old gmx_mtop_t_to_t_topology API!")
        endif()
    else()
        message("-- Gromacs not found! Gromacs-related features are disabled!")
    endif()
endif()



if(${OPENBABEL2_FOUND})
    set(BabelFiles
        babel_wrapper.h
        babel_wrapper.cpp
        pdbqt_file.h
        pdbqt_file.cpp)
endif()

add_library(pteros_io STATIC
    ${PROJECT_SOURCE_DIR}/include/pteros/core/mol_file.h
    mol_file.cpp
    vmd_molfile_plugin_wrapper.h
    vmd_molfile_plugin_wrapper.cpp

    ${BabelFiles}

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

target_link_libraries(pteros_io molfile_plugins xdrfile tng_io)
