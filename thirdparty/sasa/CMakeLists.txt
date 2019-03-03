#------------------------------
# SASA library from MDTraj
#------------------------------

add_library(sasa STATIC
    msvccompat.h
    sasa.cpp
    sasa.h    
    vectorize_sse.h
)

target_compile_options(sasa PRIVATE "-msse4.1")
target_include_directories(sasa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
