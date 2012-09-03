IF(VMD_PDB)
     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVMD_PDB")
ENDIF()

link_directories(${PROJECT_SOURCE_DIR}/lib)

# Headers of xdrfile are stored together with sources
include_directories(${PROJECT_SOURCE_DIR}/src/xdrfile)
# The same is for molfile plugins
include_directories(${PROJECT_SOURCE_DIR}/src/molfile_plugins)


add_library(pteros 
system.cpp
selection.cpp  
gnm.cpp 
grid_search.cpp 
pdb_cryst.cpp
format_recognition.cpp

${PROJECT_SOURCE_DIR}/include/pteros/core/atom.h    

${PROJECT_SOURCE_DIR}/include/pteros/core/force_field.h
#force_field.cpp

${PROJECT_SOURCE_DIR}/include/pteros/core/pteros_error.h  
${PROJECT_SOURCE_DIR}/include/pteros/core/selection.h          
${PROJECT_SOURCE_DIR}/include/pteros/core/system.h
${PROJECT_SOURCE_DIR}/include/pteros/core/grid_search.h  
${PROJECT_SOURCE_DIR}/include/pteros/core/gnm.h  

${PROJECT_SOURCE_DIR}/include/pteros/core/pdb_cryst.h
${PROJECT_SOURCE_DIR}/include/pteros/core/format_recognition.h
${PROJECT_SOURCE_DIR}/include/pteros/core/selection_parser.h
selection_parser.cpp

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

tpr_file.h
tpr_file.cpp

)

target_link_libraries(pteros molfile_plugins xdrfile)

set_target_properties(pteros PROPERTIES 
ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

