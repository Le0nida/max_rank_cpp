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
# Libspatialindex
# ----------------------------------
set(SPATIAL_INCLUDE_DIR "C:/Program Files (x86)/spatialindex/include")
set(SPATIAL_LIB_DIR "C:/Program Files (x86)/spatialindex/lib")
set(SPATIAL_BIN_DIR "C:/Program Files (x86)/spatialindex/bin")

include_directories(${SPATIAL_INCLUDE_DIR})
link_directories(${SPATIAL_LIB_DIR})

find_library(SPATIALINDEX_LIBRARY NAMES spatialindex
        PATHS "${SPATIAL_LIB_DIR}" NO_DEFAULT_PATH)
find_library(SPATIALINDEX_C_LIBRARY NAMES spatialindex_c
        PATHS "${SPATIAL_LIB_DIR}" NO_DEFAULT_PATH)

if (NOT SPATIALINDEX_LIBRARY)
    message(FATAL_ERROR "SPATIALINDEX_LIBRARY not found! Check installation paths.")
endif()
if (NOT SPATIALINDEX_C_LIBRARY)
    message(FATAL_ERROR "SPATIALINDEX_C_LIBRARY not found! Check installation paths.")
endif()

# ----------------------------------
# BOOST (for Geometry / R-Tree)
# ----------------------------------
# If you do not want to rely on environment variables, set it manually:
# set(BOOST_ROOT "C:/Libraries/boost_1_82_0")

# We can look for (at least) version 1.70, no particular COMPONENTS needed for geometry
set(Boost_INCLUDE_DIRS "C:/Program Files (x86)/C-Libraries/boost_1_87_0/boost_1_87_0")

# Add Boost include dirs
include_directories(${Boost_INCLUDE_DIRS})

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
        ${SPATIALINDEX_LIBRARY}
        ${SPATIALINDEX_C_LIBRARY}
        # Boost.Geometry is header-only, so no library needed, but if you used other boost libs,
        # you'd do something like: ${Boost_LIBRARIES}
)

# ----------------------------------
# Copy .dll files for libspatialindex
# ----------------------------------
add_custom_command(TARGET MaxRankProject POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${SPATIAL_BIN_DIR}/libspatialindex.dll"
        "${SPATIAL_BIN_DIR}/libspatialindex_c.dll"
        $<TARGET_FILE_DIR:MaxRankProject>
)
