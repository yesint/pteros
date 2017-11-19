#------------------------------
# voro++ library 
#------------------------------

add_library(voro++ SHARED
    cell.cc
    common.cc
    container.cc
    unitcell.cc
    v_compute.cc
    c_loops.cc
    v_base.cc
    wall.cc
    pre_container.cc
    container_prd.cc
)

install(TARGETS voro++
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
