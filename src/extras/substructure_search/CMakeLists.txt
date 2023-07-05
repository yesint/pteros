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

    install(TARGETS pteros_substructure_search
        LIBRARY DESTINATION lib/pteros
        ARCHIVE DESTINATION lib/pteros
    )

endif()
