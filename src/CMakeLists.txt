cmake_minimum_required(VERSION 3.6)
project(stlsvg
        VERSION 1.0
        LANGUAGES CXX)

find_package(CGAL REQUIRED COMPONENTS Core
        CMAKE_FIND_ROOT_PATH_BOTH)

add_executable(stlsvg
        stlsvg.cpp
        clipper/clipper.cpp)

set_property(TARGET stlsvg PROPERTY CXX_STANDARD 17)
target_link_libraries(stlsvg PRIVATE CGAL::CGAL)

set(CMAKE_C_OUTPUT_EXTENSION   ".o")
set(CMAKE_CXX_OUTPUT_EXTENSION ".o")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g4")
set(CMAKE_EXECUTABLE_SUFFIX    ".html")
set(CMAKE_SIZEOF_VOID_P 4)

set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(SHELL ${CMAKE_CURRENT_SOURCE_DIR}/emshell.html)
# Debugging:
#   --emrun --emit-symbol-map -s DEMANGLE_SUPPORT=1 -s SAFE_HEAP=1
set_target_properties(stlsvg PROPERTIES
        LINK_FLAGS
        "--bind --shell-file ${SHELL} --no-heap-copy -O3 -s ALLOW_MEMORY_GROWTH=1 -s MINIFY_HTML=0 -s WASM=1 -s ASYNCIFY=1")
