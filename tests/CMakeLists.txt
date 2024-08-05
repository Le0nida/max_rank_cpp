cmake_minimum_required(VERSION 3.10)

project(MaxRankProjectTests)

# Enable testing
enable_testing()

# Add GoogleTest
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/release-1.10.0.zip
        DOWNLOAD_EXTRACT_TIMESTAMP true
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Add test executable for QNode
add_executable(test_qnode test_qnode.cpp)
target_link_libraries(test_qnode gtest_main qtree_lib)
include(GoogleTest)
gtest_discover_tests(test_qnode)

# Add test executable for QTree
add_executable(test_qtree test_qtree.cpp)
target_link_libraries(test_qtree gtest_main qtree_lib)
include(GoogleTest)
gtest_discover_tests(test_qtree)
