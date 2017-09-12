#------------------------------
# xdrfile library from Gromacs
#------------------------------

add_library(xdrfile STATIC
    xdrfile.c
    xdrfile_trr.c
    xdrfile_xtc.c
    xdrfile.h
    xdrfile_trr.h
    xdrfile_xtc.h
)

