cmake_minimum_required(VERSION 3.10)
project(MaxRankProject)

set(CMAKE_CXX_STANDARD 17)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

add_executable(MaxRankProject main.cpp)
target_include_directories(MaxRankProject PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(MaxRankProject PRIVATE qtree_lib)
