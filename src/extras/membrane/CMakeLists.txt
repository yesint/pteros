add_library(pteros_membrane SHARED
    lipid_membrane.cpp
    lipid_tail_descr.cpp
    lipid_tail.cpp
    lipid_molecule.cpp
    lipid_species.cpp
    quad_surface.cpp
    local_patch.cpp
    per_species_properties.cpp
    lipid_group.cpp
    )

if(MINGW)
    # Workaround CMake problem on mingw with incorrect implicit linking
    target_link_libraries(pteros_membrane PRIVATE spdlog::spdlog)
endif()

target_link_libraries(pteros_membrane PRIVATE pteros voro2d voro3d)

set_target_properties(pteros_membrane PROPERTIES TO_INSTALL ON)



