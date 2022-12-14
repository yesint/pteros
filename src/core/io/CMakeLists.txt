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
    gmx_box_utils.h
    gmx_box_utils.cpp
)

if(WITH_TNG)
    target_sources(pteros_io PRIVATE
        tng_file.h
        tng_file.cpp)
    target_link_libraries(pteros_io PRIVATE tng_interface)
endif()

if(WITH_OPENBABEL)
    target_sources(pteros_io PRIVATE
        babel_wrapper.h
        babel_wrapper.cpp
        pdbqt_file.h
        pdbqt_file.cpp
        mol2_file.h
        mol2_file.cpp)
    target_link_libraries(pteros_io PRIVATE
        openbabel_interface
        babel_utils)
endif()

if(WITH_GROMACS)
    target_sources(pteros_io PRIVATE
        tpr_file.h
        tpr_file.cpp)
    target_link_libraries(pteros_io PRIVATE gromacs_interface)
endif()

# Unconditional libs
target_link_libraries(pteros_io PRIVATE
    pteros_interface
    molfile_plugins
    xdrfile
    spdlog::spdlog
    Eigen3::Eigen
    fmt::fmt
)
