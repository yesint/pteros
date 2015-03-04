#------------------
# pteros_rmsf
#------------------

#add_executable(pteros_rmsf
#    rmsf_program.cpp
#    rmsf.h
#    rmsf.cpp
#)
#target_link_libraries(pteros_rmsf pteros pteros_analysis ${Boost_LIBRARIES})

#------------------------
# pteros_contacts_finder
#------------------------

add_executable(pteros_contacts_finder
    contacts_finder_program.cpp
    contacts_finder.h
    contacts_finder.cpp
)
target_link_libraries(pteros_contacts_finder pteros pteros_analysis ${Boost_LIBRARIES})

#------------------------
# pteros_solvate
#------------------------

add_executable(pteros_solvate
    solvate.cpp
)
target_link_libraries(pteros_solvate pteros pteros_analysis ${Boost_LIBRARIES})


#--------------
# Installation
#--------------

install(TARGETS
    #pteros_rmsf
    pteros_contacts_finder
    pteros_solvate

    RUNTIME DESTINATION bin/tools
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
