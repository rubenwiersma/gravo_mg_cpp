cmake_minimum_required(VERSION 3.16)
project(gravomg)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Add your project files
SET(SRC_FILES
    src/multigrid_solver.cpp
    src/utility.cpp
)

SET(HEADER_FILES
    include/gravomg/multigrid_solver.h
    include/gravomg/utility.h
)

add_library(gravomg ${SRC_FILES} ${HEADER_FILES})
target_link_libraries(gravomg PUBLIC igl::core)
target_include_directories(gravomg PUBLIC include)