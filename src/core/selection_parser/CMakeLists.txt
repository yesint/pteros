#--------------------------------
# pteros selection parser library
#--------------------------------

target_sources(pteros PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/selection_macro.h
    ${CMAKE_CURRENT_LIST_DIR}/selection_parser.h
    ${CMAKE_CURRENT_LIST_DIR}/peglib.h
    ${CMAKE_CURRENT_LIST_DIR}/selection_parser.cpp
)

target_include_directories(pteros PRIVATE ${CMAKE_CURRENT_LIST_DIR})
