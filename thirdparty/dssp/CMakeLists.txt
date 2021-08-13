#---------------------
# DSSP library
#---------------------

add_library(dssp STATIC
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

target_link_libraries(dssp PRIVATE Eigen3::Eigen spdlog::spdlog fmt::fmt)
target_include_directories(dssp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
