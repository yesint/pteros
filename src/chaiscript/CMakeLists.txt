if(CMAKE_HOST_UNIX)
  SET(DYNAMIC_LOADER "dl")
endif(CMAKE_HOST_UNIX)

FIND_LIBRARY(READLINE_LIBRARY NAMES readline PATH /usr/lib /usr/local/lib /opt/local/lib)
MESSAGE(STATUS "Detecting readline support")
if (READLINE_LIBRARY)
    MESSAGE(STATUS "Found: ${READLINE_LIBRARY}")
    SET (READLINE_LIB readline)
    ADD_DEFINITIONS(/DREADLINE_AVAILABLE)
else(READLINE_LIBRARY)
    MESSAGE(STATUS "Not Found")
    SET (READLINE_LIB )
    SET (READLINE_FLAG )
endif(READLINE_LIBRARY)

include_directories(${PROJECT_SOURCE_DIR}/include/chaiscript)

add_executable(pteros_script
main.cpp
engine.cpp
${PROJECT_SOURCE_DIR}/include/pteros/bindings/chaiscript_engine.h
wrap_eigen.cpp wrap_eigen.h
wrap_system.cpp wrap_system.h
wrap_selection.cpp wrap_selection.h
)

target_link_libraries(pteros_script ${DYNAMIC_LOADER} ${READLINE_LIB} pteros pteros_sim xdrfile)
