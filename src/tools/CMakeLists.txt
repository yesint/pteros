#------------------------
# pteros_solvate
#------------------------

add_executable(pteros_solvate_program
    solvate_program.cpp
)
target_link_libraries(pteros_solvate_program PRIVATE pteros pteros_analysis pteros_solvate ${Boost_LIBRARIES})
target_include_directories(pteros_solvate_program PRIVATE pteros_solvate)

#------------------------
# pteros_rearrange
#------------------------

add_executable(pteros_rearrange_program
    rearrange.cpp
)
target_link_libraries(pteros_rearrange_program PRIVATE pteros pteros_analysis ${Boost_LIBRARIES})

#--------------
# Installation
#--------------

install(TARGETS
    pteros_solvate_program
    pteros_rearrange_program

    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
