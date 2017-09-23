cmake_minimum_required(VERSION 2.8)
project(main)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17")

set(SOURCE_FILES
main.cpp
)

include_directories(.)

add_executable(main ${SOURCE_FILES})

LINK_DIRECTORIES(
/usr/lib/x86_64-linux-gnu/
/usr/include/libpng16
)
TARGET_LINK_LIBRARIES(main pthread png16)

