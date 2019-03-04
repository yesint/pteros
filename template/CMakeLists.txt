###############################################
## This is general CMakeLists.txt file       ##
## which is needed to link your own program  ##
## with Pteros libraries                     ##
###############################################

cmake_minimum_required(VERSION 3.1)

project(pteros_template_program)

find_package(pteros REQUIRED)

add_executable(pteros_template_executable  main.cpp)

target_link_libraries(pteros_template_executable PRIVATE pteros::pteros)