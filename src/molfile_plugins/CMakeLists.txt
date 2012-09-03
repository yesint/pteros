
add_library(molfile_plugins 
molfile_plugin.h
periodic_table.h
vmdplugin.h
largefiles.h
endianswap.h
libmolfile_plugin.h

readpdb.h
pdbplugin.c

fastio.h
dcdplugin.c
)

set_target_properties(molfile_plugins PROPERTIES
ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_executable(test_plugins test.cpp)
target_link_libraries(test_plugins molfile_plugins pteros)
