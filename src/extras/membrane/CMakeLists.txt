add_library(pteros_membrane SHARED
    membrane.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/extras/membrane.h
    )

target_link_libraries(pteros_membrane spdlog::spdlog pteros voro++)

install(TARGETS pteros_membrane
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)


