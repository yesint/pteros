#add_executable(threads_test test_threads.cpp)
#target_link_libraries(threads_test
#    pteros pteros_analysis ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

add_executable(pteros_test
test_selection.cpp
    #${PROJECT_SOURCE_DIR}/src/core/peg_parser.cpp
    #${PROJECT_SOURCE_DIR}/include/pteros/core/peg_parser.h
)

#set_target_properties(pteros_test PROPERTIES CXX_STANDARD 17)
#set_property(TARGET pteros_test
#    PROPERTY CXX_STANDARD 17
#    PROPERTY CXX_STANDARD_REQUIRED ON
#)


target_link_libraries(pteros_test pteros_analysis pteros pteros_voronoi_packing)

install(TARGETS
    pteros_test

    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

#add_executable(pteros_peglib_test
#grammar2.cpp
#peglib.h
#)
#target_link_libraries(pteros_peglib_test pteros)


#add_subdirectory(${PROJECT_SOURCE_DIR}/src/python/pybind11_src)
#INCLUDE_DIRECTORIES(pybind11_src/include)
#pybind11_add_module(module1 test_1.cpp)

#add_executable(test_sim test_sim.cpp)
#target_link_libraries(test_sim  pteros)

#Installation
#install(TARGETS pteros_test
#    RUNTIME DESTINATION bin
#)

#install(TARGETS threads_test
#    RUNTIME DESTINATION bin
#)
