add_library(qtree_lib qtree.cpp geom.cpp qnode.cpp halfspace.cpp query.cpp cell.cpp maxrank.cpp main.cpp utils.cpp)
target_include_directories(qtree_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
