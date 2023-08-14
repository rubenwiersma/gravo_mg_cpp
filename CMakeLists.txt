cmake_minimum_required(VERSION 3.16)
project(gravomg)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Eigen
if(TARGET Eigen3::Eigen)
    if(${GC_EIGEN_LOCATION})
        # Change this to a FATAL_ERROR if the constraint should be observed
        message(WARNING "Target Eigen3::Eigen already exists. ")
    endif()
else()
    include(FetchContent)

    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG         3.3.7
    )

    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()

    add_library (eigen INTERFACE)
    add_library (Eigen3::Eigen ALIAS eigen)
    target_include_directories(eigen INTERFACE ${eigen_SOURCE_DIR})
endif()

# ===  done finding Eigen

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
target_link_libraries(gravomg PUBLIC Eigen3::Eigen)
target_include_directories(gravomg PUBLIC include)

# Conditionally add the -fPIC flag for shared library on supported platforms
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "amd64")
    target_compile_options(gravomg PRIVATE $<$<CXX_COMPILER_ID:GNU>:-fPIC>)
endif()
