#------------------------------
# xdrfile library from Gromacs
#------------------------------

add_library(xdrfile STATIC
    ms_stdint.h
    trr_header.h
    xdrfile.c
    xdrfile.h
    xdrfile_trr.c
    xdrfile_trr.h
    xdrfile_xtc.c
    xdrfile_xtc.h
    xdr_seek.c
    xdr_seek.h
    xdr_utils.cpp
    xdr_utils.h
)

target_include_directories(xdrfile PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
