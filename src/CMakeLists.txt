cmake_minimum_required(VERSION 3.15)
project(stlsvg
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL REQUIRED COMPONENTS Core)

add_library(clipper OBJECT
        clipper/clipper.cpp
        clipper/clipper.hpp)
target_include_directories(clipper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/clipper)

if(EMSCRIPTEN)
    # Generate index.html (and implicitly index.js, index.wasm)
    add_executable(index stlsvg.cpp)
    target_link_libraries(index PRIVATE clipper CGAL::CGAL)

    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    
    # Define shell file path
    set(SHELL_FILE ${CMAKE_CURRENT_SOURCE_DIR}/emshell.html)

    target_link_options(index PRIVATE
        "SHELL:--bind"
        "SHELL:--shell-file ${SHELL_FILE}"
        "-O3"
        "SHELL:-sALLOW_MEMORY_GROWTH=1"
        "SHELL:-sMINIFY_HTML=0"
        "SHELL:-sWASM=1"
        "SHELL:-sASYNCIFY=1"
    )
else()
    # Native testing executable
    add_executable(stlsvg-test main.cpp stlsvg.cpp)
    target_link_libraries(stlsvg-test PRIVATE clipper CGAL::CGAL)
endif()
