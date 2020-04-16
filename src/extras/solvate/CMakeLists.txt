add_library(pteros_solvate SHARED
        solvate.cpp
        ${PROJECT_SOURCE_DIR}/include/pteros/extras/solvate.h
        )

if(MINGW)
    target_link_libraries(pteros_solvate PRIVATE spdlog::spdlog)
endif()
target_link_libraries(pteros_solvate PRIVATE pteros)

install(TARGETS pteros_solvate
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
