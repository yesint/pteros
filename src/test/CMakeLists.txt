# Link with pteros libs
link_directories(${PROJECT_SOURCE_DIR}/lib)

#add_executable(threads_test test_threads.cpp)
#target_link_libraries(threads_test
#    pteros pteros_analysis ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})


add_executable(pteros_test test_selection.cpp)
target_link_libraries(pteros_test
    pteros pteros_analysis ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# Only need parent pteros include dir
include_directories(${PROJECT_SOURCE_DIR}/include)

#add_executable(test_sim test_sim.cpp)
#target_link_libraries(test_sim  pteros)

#Installation
install(TARGETS pteros_test
    RUNTIME DESTINATION bin/test
)

#install(TARGETS threads_test
#    RUNTIME DESTINATION bin/test
#)
