add_library(pteros_solvate SHARED
        solvate.cpp
        ${PROJECT_SOURCE_DIR}/include/pteros/extras/solvate.h
        )

if(MINGW)
    target_link_libraries(pteros_solvate PRIVATE spdlog::spdlog)
endif()
target_link_libraries(pteros_solvate PRIVATE pteros)

set_target_properties(pteros_solvate PROPERTIES TO_INSTALL ON)
