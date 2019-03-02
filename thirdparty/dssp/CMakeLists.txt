#---------------------
# DSSP library
#---------------------

add_library(dssp STATIC
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

target_link_libraries(dssp Eigen3::Eigen)
