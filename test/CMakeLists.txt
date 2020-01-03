cmake_minimum_required(VERSION 2.6)
message(STATUS " ***** TEST ***** ")
# We need thread support
find_package(Threads REQUIRED)

# Locate GTest
find_package(PkgConfig)
pkg_check_modules(GTEST REQUIRED gtest>=1.7.0)
pkg_check_modules(GMOCK REQUIRED gmock>=1.7.0)

include_directories(
    ${GTEST_INCLUDE_DIRS}
    ${GMOCK_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/external/include
    ${PROJECT_SOURCE_DIR}/src
)


link_directories(${PROJECT_SOURCE_DIR}/external/archive ${PROJECT_SOURCE_DIR}/external/plugin) 
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17")
message(STATUS "Using CXX flags: " ${CMAKE_CXX_FLAGS})

# Link runTests with what we want to test and the GTest and pthread library
file(GLOB SOURCE_FILES *.cpp */*.cpp)
add_executable(runTests ${SOURCE_FILES})
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread module)

