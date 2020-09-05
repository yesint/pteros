add_subdirectory(gnm)
add_subdirectory(membrane)
add_subdirectory(solvate)
# Only compile substructure search if babel is found
if(OPENBABEL2_FOUND OR OPENBABEL3_FOUND)
    add_subdirectory(substructure_search)
endif()

