#------------------
# Analysis library
#------------------

add_library(pteros_analysis SHARED    
    message_channel.h    
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/options.h
    options.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/jump_remover.h
    jump_remover.cpp    

    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/trajectory_reader.h
    trajectory_reader.cpp
    traj_file_reader.h
    traj_file_reader.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/task_base.h
    task_base.cpp
    task_driver.h
    task_driver.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/analysis/task_plugin.h
    task_plugin.cpp
    data_container.h
    )

if(MINGW)
    # Workaround CMake problem on mingw with incorrect implicit linking
    target_link_libraries(pteros_analysis PRIVATE spdlog::spdlog)
endif()

target_link_libraries(pteros_analysis PUBLIC pteros)

#--------------
# Installation
#--------------

install(TARGETS pteros_analysis EXPORT pterosTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
