add_library(pteros_voronoi_packing SHARED
    voronoi_packing.cpp
    ${PROJECT_SOURCE_DIR}/include/pteros/extras/voronoi_packing.h
    )

if(MINGW)
    # Workaround CMake problem on mingw with incorrect implicit linking
    target_link_libraries(pteros_voronoi_packing PRIVATE spdlog::spdlog)
endif()

target_link_libraries(pteros_voronoi_packing PRIVATE pteros voro++)

install(TARGETS pteros_voronoi_packing
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
