# - Config file for pteros

include(CMakeFindDependencyMacro)
find_dependency(Eigen3)
find_dependency(spdlog)
include(${CMAKE_CURRENT_LIST_DIR}/pterosTargets.cmake)

# Provide information in variables as well for backwards compatibility
set(pteros_LIBRARIES pteros::pteros)
get_target_property(pteros_INCLUDE_DIRS pteros::pteros INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(pteros_DEFINITIONS pteros::pteros INTERFACE_COMPILE_DEFINITIONS)
