find_package(PNG REQUIRED)
find_package(JPEG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
add_executable(splat
    itwom3.0.c
    splat.cpp)
target_link_libraries(splat bz2 ${PNG_LIBRARIES} ${JPEG_LIBRARIES})
install(TARGETS splat DESTINATION bin)
