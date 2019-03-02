#--------------------------------
# pteros distance_search library
#--------------------------------

add_library(pteros_distance_search STATIC
    ${PROJECT_SOURCE_DIR}/include/pteros/core/distance_search.h
    distance_search.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/distance_search_within.h
    distance_search_within.cpp

    atomic_wrapper.h
    search_utils.h
    search_utils.cpp

    distance_search_base.h
    distance_search_base.cpp

    distance_search_within_base.h
    distance_search_within_base.cpp

    distance_search_within_sel.h
    distance_search_within_sel.cpp

    distance_search_contacts.h
    distance_search_contacts.cpp

    distance_search_contacts_1sel.h
    distance_search_contacts_1sel.cpp

    distance_search_contacts_2sel.h
    distance_search_contacts_2sel.cpp
)

target_link_libraries(pteros_distance_search Eigen3::Eigen)
