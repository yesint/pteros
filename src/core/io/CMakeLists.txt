#---------------------
# pteros IO library
#---------------------

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
