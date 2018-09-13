add_subdirectory(gnm)
add_subdirectory(membrane)
# Only compile substructure search if babel is found
if(${OPENBABEL2_FOUND})
    add_subdirectory(substructure_search)
endif()

