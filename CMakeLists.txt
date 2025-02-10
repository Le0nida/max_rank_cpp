cmake_minimum_required(VERSION 3.10)
project(MaxRankProject)

set(CMAKE_CXX_STANDARD 17)

# ----------------------------------
# EIGEN
# ----------------------------------
set(EIGEN3_INCLUDE_DIR "C:/Program Files (x86)/C-Libraries/eigen-3.4.0")
include_directories(${EIGEN3_INCLUDE_DIR})

# ----------------------------------
# HIGHS
# ----------------------------------
set(HIGHS_LIB_DIR "C:/Program Files (x86)/C-Libraries/HiGHS/build/bin")
set(HIGHS_INCLUDE_DIR "C:/Program Files (x86)/C-Libraries/HiGHS")
include_directories(${HIGHS_INCLUDE_DIR})
include_directories(${HIGHS_INCLUDE_DIR}/src)
include_directories(${HIGHS_INCLUDE_DIR}/src/lp_data)
include_directories(${HIGHS_INCLUDE_DIR}/src/util)
include_directories(${HIGHS_INCLUDE_DIR}/build)

# ----------------------------------
# Subdir for your sources
# ----------------------------------
add_subdirectory(src)

# ----------------------------------
# Create main executable
# ----------------------------------
add_executable(MaxRankProject src/main.cpp)
target_include_directories(MaxRankProject PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link everything
target_link_libraries(MaxRankProject PRIVATE
        qtree_lib
        ${HIGHS_LIB_DIR}/libhighs.a
)