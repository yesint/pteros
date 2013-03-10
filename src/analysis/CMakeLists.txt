set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/analysis")

if(WRAP_PYTHON)
    IF(${CMAKE_SYSTEM_NAME} MATCHES "linux")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEXTRA_PYTHON_PATH\\(\\)=\\\"${CMAKE_BINARY_DIR}/lib\\\"")
    ENDIF()
    add_executable(python_analysis
        python_analysis_program.cpp
        task_python.cpp
        task_python.h
    )
    target_link_libraries(python_analysis pteros_analysis pteros_sim pteros ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
endif()

add_library(pteros_analysis
trajectory_processor.cpp
options_parser.cpp
rmsf.cpp
contacts_finder.cpp
selection_analysis.cpp
consumer.cpp
${PROJECT_SOURCE_DIR}/include/pteros/analysis/trajectory_processor.h
${PROJECT_SOURCE_DIR}/include/pteros/analysis/options_parser.h
${PROJECT_SOURCE_DIR}/include/pteros/analysis/rmsf.h
${PROJECT_SOURCE_DIR}/include/pteros/analysis/contacts_finder.h
${PROJECT_SOURCE_DIR}/include/pteros/analysis/selection_analysis.h
${PROJECT_SOURCE_DIR}/include/pteros/analysis/consumer.h
${PROJECT_SOURCE_DIR}/include/pteros/analysis/bilayer.h

task_box.cpp task_box.h
task_center.cpp task_center.h
task_rmsd.cpp task_rmsd.h
task_interaction_energy.cpp task_interaction_energy.h
task_distance_matr.cpp task_distance_matr.h
task_align_parts.cpp task_align_parts.h
task_covar_matr.cpp task_covar_matr.h
bilayer.cpp
)

add_executable(pteros_rmsf
    rmsf_program.cpp
)
target_link_libraries(pteros_rmsf pteros_analysis pteros_sim pteros ${Boost_LIBRARIES})

add_executable(pteros_contacts_finder
    contacts_finder_program.cpp
)
target_link_libraries(pteros_contacts_finder pteros_analysis pteros_sim pteros ${Boost_LIBRARIES})

add_executable(pteros_selection_analysis
    selection_analysis_program.cpp
)
target_link_libraries(pteros_selection_analysis pteros_analysis pteros_sim pteros ${Boost_LIBRARIES})
