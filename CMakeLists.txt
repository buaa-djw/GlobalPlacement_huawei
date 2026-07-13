cmake_minimum_required(VERSION 3.16)
project(GlobalPlacementHuawei LANGUAGES CXX)

# =========================
# C++ standard
# =========================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 你的 Limbo 当前是 std::__cxx11 ABI，因此统一使用 ABI=1
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=1)

# =========================
# Hard-coded Limbo paths
# =========================
set(LIMBO_ROOT
    "/home/djw/Desktop/GlobalPlacement_huawei/thrid_party/Limbo"
)

set(LIMBO_INCLUDE_DIR
    "${LIMBO_ROOT}"
)

set(LIMBO_BOOKSHELF_HEADER
    "${LIMBO_ROOT}/limbo/parsers/bookshelf/bison/BookshelfDriver.h"
)

set(LIMBO_BOOKSHELF_PARSER_LIBRARY
    "${LIMBO_ROOT}/build/limbo/parsers/bookshelf/bison/libbookshelfparser.a"
)

set(LIMBO_GZSTREAM_LIBRARY
    "${LIMBO_ROOT}/build/limbo/thirdparty/gzstream/libgzstream.a"
)

# =========================
# Required system libraries
# =========================
find_package(ZLIB REQUIRED)

# =========================
# Check Limbo files
# =========================
if(NOT EXISTS "${LIMBO_BOOKSHELF_HEADER}")
    message(FATAL_ERROR
        "Cannot find Limbo BookshelfDriver.h:\n"
        "  ${LIMBO_BOOKSHELF_HEADER}\n"
        "Please check LIMBO_ROOT."
    )
endif()

if(NOT EXISTS "${LIMBO_BOOKSHELF_PARSER_LIBRARY}")
    message(FATAL_ERROR
        "Cannot find libbookshelfparser.a:\n"
        "  ${LIMBO_BOOKSHELF_PARSER_LIBRARY}\n"
        "Please build Limbo Bookshelf parser first."
    )
endif()

if(NOT EXISTS "${LIMBO_GZSTREAM_LIBRARY}")
    message(FATAL_ERROR
        "Cannot find libgzstream.a:\n"
        "  ${LIMBO_GZSTREAM_LIBRARY}\n"
        "Please build Limbo gzstream first."
    )
endif()

message(STATUS "Using Limbo root: ${LIMBO_ROOT}")
message(STATUS "Using Limbo include directory: ${LIMBO_INCLUDE_DIR}")
message(STATUS "Using Limbo Bookshelf parser library: ${LIMBO_BOOKSHELF_PARSER_LIBRARY}")
message(STATUS "Using Limbo gzstream library: ${LIMBO_GZSTREAM_LIBRARY}")
message(STATUS "Using ZLIB: ${ZLIB_LIBRARIES}")

# =========================
# placement_db library
# =========================
add_library(placement_db STATIC
    src/db/PlacementDB.cpp
    src/evaluator/HPWLEvaluator.cpp
    src/parser/LimboBookshelfAdapter.cpp
)

target_include_directories(placement_db PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${LIMBO_INCLUDE_DIR}
)

target_compile_definitions(placement_db PUBLIC
    HAVE_LIMBO_BOOKSHELF
)

# 注意：
# libbookshelfparser.a 依赖 libgzstream.a
# libgzstream.a 依赖 zlib
# 链接顺序必须是：
#   libbookshelfparser.a -> libgzstream.a -> ZLIB::ZLIB
target_link_libraries(placement_db PUBLIC
    "${LIMBO_BOOKSHELF_PARSER_LIBRARY}"
    "${LIMBO_GZSTREAM_LIBRARY}"
    ZLIB::ZLIB
)

# =========================
# main executable
# =========================
add_executable(global_placer
    main.cpp
)

target_link_libraries(global_placer PRIVATE
    placement_db
    "${LIMBO_BOOKSHELF_PARSER_LIBRARY}"
    "${LIMBO_GZSTREAM_LIBRARY}"
    ZLIB::ZLIB
)

# =========================
# tests
# =========================
include(CTest)

if(BUILD_TESTING)
    add_executable(test_placement_db
        tests/test_placement_db.cpp
    )

    target_link_libraries(test_placement_db PRIVATE
        placement_db
        "${LIMBO_BOOKSHELF_PARSER_LIBRARY}"
        "${LIMBO_GZSTREAM_LIBRARY}"
        ZLIB::ZLIB
    )

    add_test(NAME placement_db_unit COMMAND test_placement_db)

    add_executable(test_hpwl
        tests/test_hpwl.cpp
    )

    target_link_libraries(test_hpwl PRIVATE
        placement_db
        "${LIMBO_BOOKSHELF_PARSER_LIBRARY}"
        "${LIMBO_GZSTREAM_LIBRARY}"
        ZLIB::ZLIB
    )

    add_test(NAME hpwl_unit COMMAND test_hpwl)
endif()
