link_libraries( pteros_analysis pteros)

#add_executable(threads_test test_threads.cpp)
#target_link_libraries(threads_test
#    pteros pteros_analysis ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

add_executable(pteros_test
test_selection.cpp
    #${PROJECT_SOURCE_DIR}/src/core/peg_parser.cpp
    #${PROJECT_SOURCE_DIR}/include/pteros/core/peg_parser.h
)
target_link_libraries(pteros_test ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${BOOST_PYTHON_LIB})

#add_executable(test_sim test_sim.cpp)
#target_link_libraries(test_sim  pteros)

#Installation
#install(TARGETS pteros_test
#    RUNTIME DESTINATION bin
#)

#install(TARGETS threads_test
#    RUNTIME DESTINATION bin
#)
