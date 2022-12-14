if(TARGET openbabel_interface)

    add_library(babel_utils STATIC
        selection_to_obmol.h
        selection_to_obmol.cpp
        )

    target_link_libraries(babel_utils
        PUBLIC
            pteros_interface
            openbabel_interface
            Eigen3::Eigen
    )

    target_include_directories(babel_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

endif()
