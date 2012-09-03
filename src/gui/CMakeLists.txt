# Find wx
SET(wxWidgets_USE_STATIC OFF)
SET(wxWidgets_USE_UNICODE ON)
FIND_PACKAGE(wxWidgets REQUIRED base core gl)
if(wxWidgets_FOUND)
  message("Building with wx")
  INCLUDE(${wxWidgets_USE_FILE})
endif()

find_package(Boost 1.37.0 REQUIRED thread signals)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  link_libraries(${Boost_LIBRARIES})
endif()

link_directories(${PROJECT_SOURCE_DIR}/lib)
link_libraries(pteros_visualizer pteros xdrfile)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(pteros_gui 
app.cpp app.h
main_window.cpp main_window.h
main_window_base.cpp main_window_base.h
)

if(wxWidgets_FOUND)
  target_link_libraries(pteros_gui ${wxWidgets_LIBRARIES} )
endif()

