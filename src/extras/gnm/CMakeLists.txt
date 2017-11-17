add_library(pteros_gnm SHARED
    gnm.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/extras/gnm.h
    )

target_link_libraries(pteros_gnm pteros)

install(TARGETS pteros_gnm
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)


