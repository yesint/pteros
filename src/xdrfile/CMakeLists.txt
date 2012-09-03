add_library(xdrfile 
xdrfile.c  xdrfile_trr.c  xdrfile_xtc.c
xdrfile.h
xdrfile_trr.h
xdrfile_xtc.h)

set_target_properties(xdrfile PROPERTIES 
ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)


