IF(VMD_PDB)
     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVMD_PDB")
ENDIF()

if(POWERSASA)
    # Set definition for conditional compile
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_POWERSASA")
    # Make list of SASA files
    SET(SASA_FILES "${PROJECT_SOURCE_DIR}/src/core/sasa/power_diagram.h" "${PROJECT_SOURCE_DIR}/src/core/sasa/power_sasa.h")
    message("POWERSASA code is used! \n This implies severe license restrictions! \n Read the comment in src/core/sasa/power_sasa.h carefully. \n If you don't agree with this license please use flag -DPOWERSASA=\"OFF\" !")
    # set includes for SASA
    include_directories(${PROJECT_SOURCE_DIR}/src/core/sasa)
    # POWERSASA requires c++11
    set_source_files_properties(${SASA_FILES} PROPERTIES COMPILE_FLAGS -std=c++0x)
else()
    SET(SASA_FILES "")
endif()

# Make list of DSSP files
SET(DSSP_FILES
"${PROJECT_SOURCE_DIR}/src/core/dssp/pteros_dssp_wrapper.cpp"
"${PROJECT_SOURCE_DIR}/src/core/dssp/pteros_dssp_wrapper.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/buffer.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/align-2d.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/matrix.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/dssp.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/mas.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/primitives-3d.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/structure.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/utils.h"
"${PROJECT_SOURCE_DIR}/src/core/dssp/dssp.cpp"
"${PROJECT_SOURCE_DIR}/src/core/dssp/mas.cpp"
"${PROJECT_SOURCE_DIR}/src/core/dssp/primitives-3d.cpp"
"${PROJECT_SOURCE_DIR}/src/core/dssp/structure.cpp"
"${PROJECT_SOURCE_DIR}/src/core/dssp/utils.cpp"
)

#DSSP requires c++11
set_source_files_properties(${DSSP_FILES} PROPERTIES COMPILE_FLAGS -std=c++0x)

# Headers of xdrfile are stored together with sources
include_directories(${PROJECT_SOURCE_DIR}/src/core/xdrfile)
# The same is for molfile plugins
include_directories(${PROJECT_SOURCE_DIR}/src/core/molfile_plugins)
# The same is for DSSP
include_directories(${PROJECT_SOURCE_DIR}/src/core/dssp)


# Main pteros library
add_library(pteros SHARED
${PROJECT_SOURCE_DIR}/include/pteros/pteros.h
${PROJECT_SOURCE_DIR}/include/pteros/core/pteros_error.h
${PROJECT_SOURCE_DIR}/include/pteros/core/atom.h
${PROJECT_SOURCE_DIR}/include/pteros/core/force_field.h
${PROJECT_SOURCE_DIR}/include/pteros/core/system.h
system.cpp
${PROJECT_SOURCE_DIR}/include/pteros/core/selection.h
selection_macro.h
selection.cpp  
${PROJECT_SOURCE_DIR}/include/pteros/core/selection_parser.h
selection_parser.cpp
${PROJECT_SOURCE_DIR}/include/pteros/core/gnm.h
gnm.cpp 
${PROJECT_SOURCE_DIR}/include/pteros/core/grid_search.h
grid_search.cpp 
${PROJECT_SOURCE_DIR}/include/pteros/core/pdb_cryst.h
pdb_cryst.cpp
${PROJECT_SOURCE_DIR}/include/pteros/core/format_recognition.h
format_recognition.cpp
# IO
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
pttop_file.h
pttop_file.cpp

# xdrfile
${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile.c
${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile_trr.c
${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile_xtc.c
${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile.h
${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile_trr.h
${PROJECT_SOURCE_DIR}/src/core/xdrfile/xdrfile_xtc.h

# molfile_plugins
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/molfile_plugin.h
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/periodic_table.h
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/vmdplugin.h
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/largefiles.h
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/endianswap.h
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/libmolfile_plugin.h
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/readpdb.h
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/pdbplugin.c
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/fastio.h
${PROJECT_SOURCE_DIR}/src/core/molfile_plugins/dcdplugin.c

#SASA (will be empty if not used)
${SASA_FILES}

#DSSP
${DSSP_FILES}
)

#Installation
install(TARGETS pteros
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
