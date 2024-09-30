cmake_minimum_required(VERSION 3.10)
project(MaxRankProject)

set(CMAKE_CXX_STANDARD 17)

# Set the path to Eigen
set(EIGEN3_INCLUDE_DIR "C:/Program Files (x86)/C-Libraries/eigen-3.4.0")

# Include Eigen
include_directories(${EIGEN3_INCLUDE_DIR})

# Include HIGHS
set(HIGHS_LIB_DIR "C:/Program Files (x86)/C-Libraries/HiGHS/build/bin")
set(HIGHS_INCLUDE_DIR "C:/Program Files (x86)/C-Libraries/HiGHS")
include_directories(${HIGHS_INCLUDE_DIR})
include_directories(${HIGHS_INCLUDE_DIR}/src)
include_directories(${HIGHS_INCLUDE_DIR}/src/lp_data)
include_directories(${HIGHS_INCLUDE_DIR}/src/util)
include_directories(${HIGHS_INCLUDE_DIR}/build)

# Include LevelDB
set(LEVELDB_INCLUDE_DIR "C:/Program Files (x86)/C-Libraries/leveldb/include")
set(LEVELDB_LIB_DIR "C:/Program Files (x86)/C-Libraries/leveldb/build")

include_directories(${LEVELDB_INCLUDE_DIR})

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

# Create executable for MaxRankProject
add_executable(MaxRankProject src/main.cpp)

# Include project directories
target_include_directories(MaxRankProject PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link HIGHS and LevelDB libraries
target_link_libraries(MaxRankProject PRIVATE qtree_lib
        ${HIGHS_LIB_DIR}/libhighs.a
        ${LEVELDB_LIB_DIR}/libleveldb.a
)
