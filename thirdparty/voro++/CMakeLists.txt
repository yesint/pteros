#------------------------------
# voro++ library 
#------------------------------

add_library(voro++ STATIC
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

target_include_directories(voro++ PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS voro++
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
