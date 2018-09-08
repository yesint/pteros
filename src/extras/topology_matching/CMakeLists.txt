add_library(pteros_topmatch SHARED
    topmatch.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/extras/topmatch.h
    )

target_link_libraries(pteros_topmatch pteros)

install(TARGETS pteros_topmatch
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)


