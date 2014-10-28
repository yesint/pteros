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

find_package(CUDA)

# nvcc flags
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3 -gencode arch=compute_20,code=sm_20)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")

#cuda_compile(CUDA_O cuda.cu)

#CUDA_ADD_EXECUTABLE(cuda_test ${CUDA_O})
#CUDA_ADD_LIBRARY(cuda_adaptor cuda_pteros_adaptor.cu SHARED)
#CUDA_ADD_EXECUTABLE(cuda_test cuda.cu)

#target_link_libraries(cuda_test cuda_adaptor)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
