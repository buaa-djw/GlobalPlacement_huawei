cmake_minimum_required(VERSION 3.16)
project(GlobalPlacementHuawei LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(DEFAULT_LIMBO_ROOT "/home/djw/Desktop/GlobalPlacement_huawei/thrid_party/Limbo/limbo")
set(LIMBO_ROOT "${DEFAULT_LIMBO_ROOT}" CACHE PATH "Path to Limbo source/include directory or installation prefix")

add_library(placement_db
    src/db/PlacementDB.cpp
    src/evaluator/HPWLEvaluator.cpp
    src/parser/LimboBookshelfAdapter.cpp
)
target_include_directories(placement_db PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(LIMBO_HEADER_FOUND FALSE)
set(LIMBO_INCLUDE_DIR "")
if(EXISTS "${LIMBO_ROOT}/include/limbo/parsers/bookshelf/bison/BookshelfDriver.h")
    set(LIMBO_HEADER_FOUND TRUE)
    set(LIMBO_INCLUDE_DIR "${LIMBO_ROOT}/include")
elseif(EXISTS "${LIMBO_ROOT}/limbo/parsers/bookshelf/bison/BookshelfDriver.h")
    set(LIMBO_HEADER_FOUND TRUE)
    set(LIMBO_INCLUDE_DIR "${LIMBO_ROOT}")
elseif(EXISTS "${LIMBO_ROOT}/parsers/bookshelf/bison/BookshelfDriver.h")
    set(LIMBO_HEADER_FOUND TRUE)
    get_filename_component(LIMBO_INCLUDE_DIR "${LIMBO_ROOT}/.." ABSOLUTE)
endif()

if(LIMBO_HEADER_FOUND)
    message(STATUS "Using Limbo include directory: ${LIMBO_INCLUDE_DIR}")
    target_include_directories(placement_db PRIVATE "${LIMBO_INCLUDE_DIR}")
    target_compile_definitions(placement_db PRIVATE HAVE_LIMBO_BOOKSHELF)

    set(_limbo_library_hints
        "${LIMBO_ROOT}/lib"
        "${LIMBO_ROOT}/../lib"
        "${LIMBO_ROOT}/build/lib"
        "${LIMBO_ROOT}/../build/lib"
        "${LIMBO_ROOT}/lib64"
        "${LIMBO_ROOT}/../lib64"
    )
    find_library(LIMBO_BOOKSHELF_PARSER_LIBRARY
        NAMES bookshelfparser libbookshelfparser
        HINTS ${_limbo_library_hints}
    )
    if(LIMBO_BOOKSHELF_PARSER_LIBRARY)
        message(STATUS "Using Limbo Bookshelf parser library: ${LIMBO_BOOKSHELF_PARSER_LIBRARY}")
        target_link_libraries(placement_db PUBLIC "${LIMBO_BOOKSHELF_PARSER_LIBRARY}")
    else()
        message(FATAL_ERROR "Found Limbo headers at ${LIMBO_INCLUDE_DIR}, but could not find libbookshelfparser. Build/install Limbo or set CMAKE_LIBRARY_PATH to the directory containing libbookshelfparser.a/.so.")
    endif()
else()
    message(WARNING "Limbo BookshelfDriver.h was not found under LIMBO_ROOT='${LIMBO_ROOT}'. The compatibility reader will be used. For djw's VM, configure with -DLIMBO_ROOT=/home/djw/Desktop/GlobalPlacement_huawei/thrid_party/Limbo/limbo")
endif()

add_executable(global_placer main.cpp)
target_link_libraries(global_placer PRIVATE placement_db)

include(CTest)
if(BUILD_TESTING)
    add_executable(test_placement_db tests/test_placement_db.cpp)
    target_link_libraries(test_placement_db PRIVATE placement_db)
    add_test(NAME placement_db_unit COMMAND test_placement_db)

    add_executable(test_hpwl tests/test_hpwl.cpp)
    target_link_libraries(test_hpwl PRIVATE placement_db)
    add_test(NAME hpwl_unit COMMAND test_hpwl)
endif()
