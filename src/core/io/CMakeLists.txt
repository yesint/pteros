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
    mol2_file.h
    mol2_file.cpp
    xyz_file.h
    xyz_file.cpp
    trr_file.h
    trr_file.cpp
    xtc_file.h
    xtc_file.cpp
)

if(WITH_TNGIO)
    target_sources(pteros_io PRIVATE
        tng_file.h
        tng_file.cpp)
    target_compile_definitions(pteros_io PRIVATE USE_TNGIO)
    target_link_libraries(pteros_io PRIVATE tng_io)
endif()

if(WITH_OPENBABEL AND (OPENBABEL2_FOUND OR OPENBABEL3_FOUND))
    target_sources(pteros_io PRIVATE
        babel_wrapper.h
        babel_wrapper.cpp
        pdbqt_file.h
        pdbqt_file.cpp)
    target_compile_definitions(pteros_io PRIVATE USE_OPENBABEL)
    target_include_directories(pteros_io PRIVATE ${OPENBABEL2_INCLUDE_DIR} ${OPENBABEL3_INCLUDE_DIR})
    target_link_libraries(pteros_io PRIVATE pteros_babel_utils)
endif()


if(WITH_GROMACS AND GROMACS_SOURCES AND GROMACS_LIBRARIES)
    target_compile_definitions(pteros_io PRIVATE USE_GROMACS)
    # In case of Gromacs 2020+ we need some headers from gmx source tree
    # since they are not installed any more, so do it also for earlier versions
    # to keep things simpler
    target_include_directories(pteros_io PRIVATE ${GROMACS_SOURCES}/src)
    target_include_directories(pteros_io PRIVATE ${GROMACS_SOURCES}/src/external/thread_mpi/include)

    target_sources(pteros_io PRIVATE
        tpr_file.h
        tpr_file.cpp)    
else()
    # Without gromacs link to xdrfile
    target_link_libraries(pteros_io PRIVATE xdrfile)
endif()

# Unconditional libs
target_link_libraries(pteros_io PRIVATE
    pteros_gromacs_utils
    molfile_plugins
    spdlog::spdlog
    Boost::boost
    Eigen3::Eigen)
