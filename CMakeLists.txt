cmake_minimum_required(VERSION 3.16)
project(GlobalPlacementHuawei LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(LIMBO_ROOT "" CACHE PATH "Path to Limbo installation")

add_library(placement_db
    src/db/PlacementDB.cpp
    src/parser/LimboBookshelfAdapter.cpp
)
target_include_directories(placement_db PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(LIMBO_ROOT)
    target_include_directories(placement_db PRIVATE ${LIMBO_ROOT})
    if(EXISTS "${LIMBO_ROOT}/limbo/parsers/bookshelf/bison/BookshelfDriver.h")
        target_compile_definitions(placement_db PRIVATE HAVE_LIMBO_BOOKSHELF)
    elseif(EXISTS "${LIMBO_ROOT}/parsers/bookshelf/bison/BookshelfDriver.h")
        target_include_directories(placement_db PRIVATE "${LIMBO_ROOT}/..")
        target_compile_definitions(placement_db PRIVATE HAVE_LIMBO_BOOKSHELF)
    else()
        message(WARNING "LIMBO_ROOT was set but BookshelfDriver.h was not found; compatibility reader will be used")
    endif()
endif()

add_executable(global_placer main.cpp)
target_link_libraries(global_placer PRIVATE placement_db)

include(CTest)
if(BUILD_TESTING)
    add_executable(test_placement_db tests/test_placement_db.cpp)
    target_link_libraries(test_placement_db PRIVATE placement_db)
    add_test(NAME placement_db_unit COMMAND test_placement_db)
endif()
