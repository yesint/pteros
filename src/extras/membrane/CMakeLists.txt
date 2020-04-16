add_library(pteros_membrane SHARED
    membrane.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/extras/membrane.h
    )

if(MINGW)
    target_link_libraries(pteros_membrane PRIVATE spdlog::spdlog)
endif()
target_link_libraries(pteros_membrane PRIVATE pteros voro++)

install(TARGETS pteros_membrane
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)


