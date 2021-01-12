#--------------------------------
# pteros distance_search library
#--------------------------------

target_sources(pteros PRIVATE
    ${PROJECT_SOURCE_DIR}/include/pteros/core/distance_search.h
    ${CMAKE_CURRENT_LIST_DIR}/distance_search.cpp

    ${PROJECT_SOURCE_DIR}/include/pteros/core/distance_search_within.h
    ${CMAKE_CURRENT_LIST_DIR}/distance_search_within.cpp

    ${CMAKE_CURRENT_LIST_DIR}/distance_search_base.h
    ${CMAKE_CURRENT_LIST_DIR}/distance_search_base.cpp

    ${CMAKE_CURRENT_LIST_DIR}/distance_search_within_base.h
    ${CMAKE_CURRENT_LIST_DIR}/distance_search_within_base.cpp

    ${CMAKE_CURRENT_LIST_DIR}/distance_search_within_sel.h
    ${CMAKE_CURRENT_LIST_DIR}/distance_search_within_sel.cpp

    ${CMAKE_CURRENT_LIST_DIR}/distance_search_contacts.h
    ${CMAKE_CURRENT_LIST_DIR}/distance_search_contacts.cpp

    ${CMAKE_CURRENT_LIST_DIR}/distance_search_contacts_1sel.h
    ${CMAKE_CURRENT_LIST_DIR}/distance_search_contacts_1sel.cpp

    ${CMAKE_CURRENT_LIST_DIR}/distance_search_contacts_2sel.h
    ${CMAKE_CURRENT_LIST_DIR}/distance_search_contacts_2sel.cpp
)
