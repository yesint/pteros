if(TARGET openbabel_interface)

    add_library(pteros_substructure_search SHARED
        substructure_search.cpp
        ${PROJECT_SOURCE_DIR}/include/pteros/extras/substructure_search.h
        )

    target_link_libraries(pteros_substructure_search
        PRIVATE
            pteros
            babel_utils
            openbabel_interface
    )

    set_target_properties(pteros_substructure_search PROPERTIES TO_INSTALL ON)
endif()
