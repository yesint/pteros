#------------------
# Analysis library
#------------------

add_library(pteros_analysis SHARED    
    message_channel.h
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/trajectory_processor.h
    trajectory_processor.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/consumer_base.h
    consumer_base.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/consumer.h
    consumer.cpp    
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/options.h
    options.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/jump_remover.h
    jump_remover.cpp    

    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/trajectory_reader.h
    trajectory_reader.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/task_base.h
    task_base.cpp
    data_container.h
    )

target_link_libraries(pteros_analysis pteros ${Boost_LIBRARIES})

#--------------
# Installation
#--------------

install(TARGETS pteros_analysis
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
