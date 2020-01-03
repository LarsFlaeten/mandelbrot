set(LIBNAME "module")

# Sources
message(STATUS "****   Building library: " ${LIBNAME} " ****")
message(STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR})
set(LIB_SOURCES "${PROJECT_SOURCE_DIR}/src" CACHE INTERNAL "")
set(LIB_VERSION_MAJOR 1)
set(LIB_VERSION 0)

include_directories(${LIB_SOURCES} ${PROJECT_SOURCE_DIR}/external/plugin ${PROJECT_SOURCE_DIR}/external/archive)

file(GLOB SOURCE_FILES *.cpp */*.cpp)

# Static or shared?
set(LIBTYPE STATIC)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib/archive)
if(BUILD_SHARED)
	set(LIBTYPE SHARED)
    set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib/plugin)
endif(BUILD_SHARED)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17")
set (LIBS "")

add_library(${LIBNAME} ${LIBTYPE} ${SOURCE_FILES})


target_link_libraries(${LIBNAME} ${LIBS})

if(BUILD_SHARED)
    # Adds SO Version and subversion. To be added to ensure ABI/API compatibility.
    SET_TARGET_PROPERTIES(${LIBNAME} PROPERTIES SOVERSION ${LIB_VERSION_MAJOR} VERSION ${LIB_VERSION})
endif(BUILD_SHARED)

link_directories(${LIBRARY_OUTPUT_PATH})

