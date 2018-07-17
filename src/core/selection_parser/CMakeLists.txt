#--------------------------------
# pteros selection parser library
#--------------------------------

add_library(pteros_selection_parser STATIC
    selection_macro.h
    selection_parser.h
    peglib.h
    selection_parser.cpp
)

# Add path to molfile plugins directory to use VMD periodic table routines
target_include_directories(pteros_selection_parser PRIVATE ${PROJECT_SOURCE_DIR}/thirdparty/molfile_plugins)
