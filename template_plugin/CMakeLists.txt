cmake_minimum_required(VERSION 3.1)

project(pteros_analysis_plugin)

find_package(pteros REQUIRED)
find_package(pybind11 REQUIRED)

###############################################
## You have to set these variables:          ##
###############################################

# Name of your plugin. This becomes the name of Python module!
SET(PLUGIN_NAME example_plugin)

# The list of source files for this plugin
SET(PLUGIN_FILES plugin_main.cpp)

###############################################
## Do not modify below!                      ##
###############################################

set_source_files_properties(${PLUGIN_FILES} PROPERTIES COMPILE_DEFINITIONS PLUGIN_NAME=${PLUGIN_NAME})

# Make python extension module
pybind11_add_module(plugin ${PLUGIN_FILES})

target_link_libraries(plugin PRIVATE pteros::pteros pteros::pteros_analysis)

set_target_properties(plugin PROPERTIES OUTPUT_NAME ${PLUGIN_NAME})

#Installation
install(TARGETS plugin
    LIBRARY DESTINATION ${pteros_PLUGINS_DIR}
    ARCHIVE DESTINATION ${pteros_PLUGINS_DIR}
)
