#------------------------
# pteros_solvate
#------------------------

add_executable(pteros_solvate_program
    solvate_program.cpp
)
if(MINGW)
    target_link_libraries(pteros_solvate_program PRIVATE spdlog::spdlog)
endif()
target_link_libraries(pteros_solvate_program PRIVATE pteros pteros_analysis pteros_solvate)
target_include_directories(pteros_solvate_program PRIVATE pteros_solvate)
set_target_properties(pteros_solvate_program PROPERTIES OUTPUT_NAME pteros_solvate)

# install
set_target_properties(pteros_solvate_program PROPERTIES TO_INSTALL ON)

#------------------------
# pteros_rearrange
#------------------------

add_executable(pteros_rearrange_program
    rearrange_program.cpp
)
if(MINGW)
    target_link_libraries(pteros_rearrange_program PRIVATE spdlog::spdlog)
endif()
target_link_libraries(pteros_rearrange_program PRIVATE pteros pteros_analysis)
set_target_properties(pteros_rearrange_program PROPERTIES OUTPUT_NAME pteros_rearrange)

# Install
set_target_properties(pteros_rearrange_program PROPERTIES TO_INSTALL ON)
