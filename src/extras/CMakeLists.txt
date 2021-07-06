add_subdirectory(gnm)
add_subdirectory(membrane)
add_subdirectory(solvate)
add_subdirectory(voronoi_packing)
# Only compile substructure search if babel is found
if(OPENBABEL2_FOUND OR OPENBABEL3_FOUND)
    add_subdirectory(substructure_search)
endif()

