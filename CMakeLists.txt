cmake_minimum_required(VERSION 3.16)
project(GlobalPlacementHuawei LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=1)

set(LIMBO_ROOT "/home/djw/Desktop/GlobalPlacement_huawei/thrid_party/Limbo" CACHE PATH "Limbo root directory")
set(LIMBO_INCLUDE_DIR "${LIMBO_ROOT}")
set(LIMBO_BOOKSHELF_HEADER "${LIMBO_ROOT}/limbo/parsers/bookshelf/bison/BookshelfDriver.h")
set(LIMBO_BOOKSHELF_PARSER_LIBRARY "${LIMBO_ROOT}/build/limbo/parsers/bookshelf/bison/libbookshelfparser.a")
set(LIMBO_GZSTREAM_LIBRARY "${LIMBO_ROOT}/build/limbo/thirdparty/gzstream/libgzstream.a")

find_package(ZLIB REQUIRED)

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
message(STATUS "Using ZLIB: ${ZLIB_LIBRARIES}")

add_library(placement_core STATIC
  src/db/PlacementDB.cpp
  src/parser/LimboBookshelfAdapter.cpp
  src/evaluator/HPWLEvaluator.cpp
  src/density/ExactBinGrid.cpp
  src/density/ExactDensityEvaluator.cpp
  src/density/PaperOverlapFunction.cpp
  src/density/ExactDensitySubgradient.cpp
  src/wirelength/ExactWirelengthSubgradient.cpp
  src/objective/NonsmoothObjectiveEvaluator.cpp
)

target_include_directories(placement_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${LIMBO_INCLUDE_DIR})
if(HAVE_LIMBO_FILES)
  target_compile_definitions(placement_core PUBLIC HAVE_LIMBO_BOOKSHELF)
  target_link_libraries(placement_core PUBLIC "${LIMBO_BOOKSHELF_PARSER_LIBRARY}" "${LIMBO_GZSTREAM_LIBRARY}" ZLIB::ZLIB)
else()
  target_link_libraries(placement_core PUBLIC ZLIB::ZLIB)
endif()

add_executable(global_placer main.cpp)
target_link_libraries(global_placer PRIVATE placement_core)

include(CTest)
if(BUILD_TESTING)
  foreach(t placement_db hpwl exact_bin_grid exact_density_evaluator exact_density_subgradient exact_wirelength_subgradient nonsmooth_objective adaptec2_smoke)
    add_executable(test_${t} tests/test_${t}.cpp)
    target_link_libraries(test_${t} PRIVATE placement_core)
    add_test(NAME ${t}_unit COMMAND test_${t})
  endforeach()
endif()
