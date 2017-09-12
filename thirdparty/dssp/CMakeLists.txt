#---------------------
# DSSP library
#---------------------

add_library(dssp SHARED
    align-2d.h
    matrix.h
    dssp.h
    mas.h
    primitives-3d.h
    structure.h
    utils.h
    dssp.cpp
    mas.cpp
    primitives-3d.cpp
    structure.cpp
    utils.cpp
)

install(TARGETS dssp
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
