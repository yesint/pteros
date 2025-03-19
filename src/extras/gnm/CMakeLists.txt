add_library(pteros_gnm SHARED
    gnm.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/extras/gnm.h
    )

if(MINGW)
    # Workaround CMake problem on mingw with incorrect implicit linking
    target_link_libraries(pteros_gnm PRIVATE spdlog::spdlog)
endif()

target_link_libraries(pteros_gnm PRIVATE pteros)

#list(APPEND PTEROS_TARGETS_TO_INSTALL pteros_gnm)
set_target_properties(pteros_gnm PROPERTIES TO_INSTALL ON)