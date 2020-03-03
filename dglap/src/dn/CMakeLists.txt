# dn because of Duff Neill - the author
set(MODULE_NAME dn)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SOURCES_LIB ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
string(REPLACE ".cxx" ".hh" HEADERS_LIB "${SOURCES_LIB}")

message(STATUS "SOURCES_LIB: ${SOURCES_LIB}")
message(STATUS "HEADERS_LIB: ${HEADERS_LIB}")

add_library(${NAME_LIB} SHARED ${SOURCES_LIB})
target_include_directories(${NAME_LIB} PUBLIC ${GSL_INCLUDE_DIRS})
target_link_libraries(${NAME_LIB} PUBLIC ${GSL_LIBRARIES})

# set(SWIG_TARGET_LINK_LIBRARIES ${GSL_LIBRARIES} ${NAME_LIB})
swigify()
install_targets_custom()
