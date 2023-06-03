#------------------------------
# voro++ library 
#------------------------------

add_library(voro3d STATIC ${CMAKE_CURRENT_SOURCE_DIR}/3d/voro++.cc)
target_include_directories(voro3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/3d)

add_library(voro2d STATIC ${CMAKE_CURRENT_SOURCE_DIR}/2d/voro++_2d.cc)
target_include_directories(voro2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/2d)
