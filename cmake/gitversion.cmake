cmake_minimum_required(VERSION 3.0.0)

message(STATUS "Resolving GIT Version")

set(_build_version "unknown")

find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY "${work_dir}"
    OUTPUT_VARIABLE _build_version
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message( STATUS "GIT hash: ${_build_version}")
else()
  message(STATUS "GIT not found, no version info available")
endif()

string(TIMESTAMP _time_stamp)

configure_file(${work_dir}/src/core/version.cpp.in ${work_dir}/src/core/version.cpp @ONLY)

