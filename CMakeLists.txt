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
    CACHE PATH "Limbo root directory"
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
set(HAVE_LIMBO_FILES TRUE)
if(NOT EXISTS "${LIMBO_BOOKSHELF_HEADER}")
    message(WARNING "Cannot find Limbo BookshelfDriver.h: ${LIMBO_BOOKSHELF_HEADER}; using compatibility parser")
    set(HAVE_LIMBO_FILES FALSE)
endif()
if(NOT EXISTS "${LIMBO_BOOKSHELF_PARSER_LIBRARY}")
    message(WARNING "Cannot find libbookshelfparser.a: ${LIMBO_BOOKSHELF_PARSER_LIBRARY}; using compatibility parser")
    set(HAVE_LIMBO_FILES FALSE)
endif()
if(NOT EXISTS "${LIMBO_GZSTREAM_LIBRARY}")
    message(WARNING "Cannot find libgzstream.a: ${LIMBO_GZSTREAM_LIBRARY}; using compatibility parser")
    set(HAVE_LIMBO_FILES FALSE)
endif()

message(STATUS "Using Limbo root: ${LIMBO_ROOT}")
message(STATUS "Using Limbo include directory: ${LIMBO_INCLUDE_DIR}")
message(STATUS "Using ZLIB: ${ZLIB_LIBRARIES}")
if(HAVE_LIMBO_FILES)
    message(STATUS "Using Limbo Bookshelf parser library: ${LIMBO_BOOKSHELF_PARSER_LIBRARY}")
    message(STATUS "Using Limbo gzstream library: ${LIMBO_GZSTREAM_LIBRARY}")
endif()

# =========================
# placement_db library
# =========================
add_library(placement_db STATIC
    src/db/PlacementDB.cpp
    src/evaluator/HPWLEvaluator.cpp
    src/evaluator/DensityEvaluator.cpp
    src/evaluator/ObjectiveEvaluator.cpp
    src/optimizer/GlobalPlacer.cpp
    src/parser/LimboBookshelfAdapter.cpp
    src/grid/BinGrid.cpp
    src/utils/Logger.cpp
)

target_include_directories(placement_db PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${LIMBO_INCLUDE_DIR}
)

if(HAVE_LIMBO_FILES)
    target_compile_definitions(placement_db PUBLIC HAVE_LIMBO_BOOKSHELF)
    target_link_libraries(placement_db PUBLIC
        "${LIMBO_BOOKSHELF_PARSER_LIBRARY}"
        "${LIMBO_GZSTREAM_LIBRARY}"
        ZLIB::ZLIB
    )
else()
    target_link_libraries(placement_db PUBLIC ZLIB::ZLIB)
endif()

# =========================
# main executable
# =========================
add_executable(global_placer
    main.cpp
)

target_link_libraries(global_placer PRIVATE
    placement_db
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
    )

    add_test(NAME placement_db_unit COMMAND test_placement_db)

    add_executable(test_hpwl
        tests/test_hpwl.cpp
    )

    target_link_libraries(test_hpwl PRIVATE
        placement_db
    )

    add_test(NAME hpwl_unit COMMAND test_hpwl)

    add_executable(test_bin_grid
        tests/test_bin_grid.cpp
    )
    target_link_libraries(test_bin_grid PRIVATE placement_db)
    add_test(NAME bin_grid_unit COMMAND test_bin_grid)

    add_executable(test_logger
        tests/test_logger.cpp
    )
    target_link_libraries(test_logger PRIVATE placement_db)
    add_test(NAME logger_unit COMMAND test_logger)

    add_executable(test_density_evaluator
        tests/test_density_evaluator.cpp
    )
    target_link_libraries(test_density_evaluator PRIVATE placement_db)
    add_test(NAME density_evaluator_unit COMMAND test_density_evaluator)

    add_executable(test_objective_evaluator tests/test_objective_evaluator.cpp)
    target_link_libraries(test_objective_evaluator PRIVATE placement_db)
    add_test(NAME objective_evaluator_unit COMMAND test_objective_evaluator)

    add_executable(test_global_placer tests/test_global_placer.cpp)
    target_link_libraries(test_global_placer PRIVATE placement_db)
    add_test(NAME global_placer_unit COMMAND test_global_placer)
endif()
