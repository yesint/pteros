#------------------------
# pteros_solvate
#------------------------

add_executable(pteros_solvate
    solvate.cpp
)
 target_link_libraries(pteros_solvate pteros pteros_analysis ${Boost_LIBRARIES})

#------------------------
# pteros_rearrange
#------------------------

add_executable(pteros_rearrange
    rearrange.cpp
)
target_link_libraries(pteros_rearrange pteros pteros_analysis ${Boost_LIBRARIES})

#--------------
# Installation
#--------------

install(TARGETS
    pteros_solvate
    pteros_rearrange

    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
