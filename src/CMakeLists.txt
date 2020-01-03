#add_subdirectory("module")

# Add defines for exe build
# add_definitions("-DXXX=")


# Sources
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/external/include" "${CMAKE_CURRENT_SOURCE_DIR}" )


set(EXENAME mandelbrot01)
file(GLOB SOURCE_FILES mandelbrot01.cpp)
add_executable(${EXENAME} ${SOURCE_FILES})
set(LIBS "mork")
target_link_libraries(${EXENAME} ${LIBS})



set(EXENAME mandelbrot02-simd-implicit)
file(GLOB SOURCE_FILES mandelbrot02-simd-implicit.cpp)
add_executable(${EXENAME} ${SOURCE_FILES})
set(LIBS "mork")
target_link_libraries(${EXENAME} ${LIBS})

set(EXENAME mandelbrot03-threads)
file(GLOB SOURCE_FILES mandelbrot03-threads.cpp)
add_executable(${EXENAME} ${SOURCE_FILES})
set(LIBS "mork")
find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
    target_compile_options(${EXENAME} PUBLIC "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(${EXENAME} ${LIBS} "${CMAKE_THREAD_LIBS_INIT}" )
else()
    target_link_libraries(${EXENAME} ${LIBS})
endif()


