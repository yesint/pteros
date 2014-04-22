#------------------
# Analysis library
#------------------

add_library(pteros_analysis SHARED
    trajectory_processor.cpp
    options_parser.cpp
    rmsf.cpp
    contacts_finder.cpp
    consumer_base.cpp
    consumer.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/message_channel.h

    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/trajectory_processor.h
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/options_parser.h
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/rmsf.h
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/contacts_finder.h
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/consumer_base.h
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/consumer.h

    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/bilayer.h
    bilayer.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/options.h
    ${PROJECT_SOURCE_DIR}/src/analysis/options.cpp
)

#------------------
# pteros_rmsf
#------------------

add_executable(pteros_rmsf
    rmsf_program.cpp
)
target_link_libraries(pteros_rmsf pteros pteros_analysis ${Boost_LIBRARIES})

#------------------------
# pteros_contacts_finder
#------------------------

add_executable(pteros_contacts_finder
    contacts_finder_program.cpp
)
target_link_libraries(pteros_contacts_finder pteros pteros_analysis ${Boost_LIBRARIES})

#--------------
# Installation
#--------------

install(TARGETS pteros_analysis pteros_rmsf pteros_contacts_finder
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
