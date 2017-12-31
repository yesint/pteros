#------------------------------
# SASA library from MDTraj
#------------------------------

add_library(sasa SHARED
    msvccompat.h
    sasa.cpp
    sasa.h    
    vectorize_sse.h
)

target_compile_options(sasa PRIVATE "-msse4.1")

install(TARGETS sasa
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
