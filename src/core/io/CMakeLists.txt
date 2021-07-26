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

if(WITH_OPENBABEL AND (OPENBABEL2_FOUND OR OPENBABEL3_FOUND))
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

if(WITH_GROMACS AND GROMACS_SOURCES AND GROMACS_LIBRARIES)
    target_compile_definitions(pteros_io PRIVATE USE_GROMACS)
    # In case of Gromacs 2020+ we need some headers from gmx source tree
    # since they are not installed any more, so do it also for earlier versions
    # to keep things simpler
    target_include_directories(pteros_io PRIVATE ${GROMACS_SOURCES}/src)
    target_sources(pteros_io PRIVATE
        tpr_file.h
        tpr_file.cpp)
endif()

# Unconditional libs
target_link_libraries(pteros_io PRIVATE
    pteros_gromacs_utils
    molfile_plugins
    xdrfile
    spdlog::spdlog
    Boost::boost
    Eigen3::Eigen
)
