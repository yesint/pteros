#--------------------------------
# pteros selection parser library
#--------------------------------

add_library(pteros_selection_parser STATIC
    selection_macro.h
    selection_parser.h
    peglib.h
    selection_parser.cpp
)

target_link_libraries(pteros_selection_parser PRIVATE pteros_core Eigen3::Eigen)
target_include_directories(pteros_selection_parser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
