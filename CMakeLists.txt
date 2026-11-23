cmake_minimum_required(VERSION 3.16)
project(pkd_discovery VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(ZLIB REQUIRED)

add_library(pkdcore
  src/util.cpp
  src/pdf_text.cpp
  src/pdf_writer.cpp
  src/knowledge_base.cpp
  src/monograph.cpp
  src/pk.cpp
  src/repurpose.cpp
  src/validate.cpp
)
target_include_directories(pkdcore PUBLIC include)
target_link_libraries(pkdcore PUBLIC ZLIB::ZLIB)
target_compile_definitions(pkdcore PUBLIC
  PKD_DATA_DIR="${CMAKE_SOURCE_DIR}/data"
  PKD_SAMPLES_DIR="${CMAKE_SOURCE_DIR}/samples")
target_compile_options(pkdcore PRIVATE -Wall -Wextra)

add_executable(pkd src/main.cpp)
target_link_libraries(pkd PRIVATE pkdcore)

add_executable(mkpdf tools/mkpdf.cpp)
target_link_libraries(mkpdf PRIVATE pkdcore)

enable_testing()
add_executable(pkd_tests tests/test_main.cpp)
target_link_libraries(pkd_tests PRIVATE pkdcore)
add_test(NAME unit_tests COMMAND pkd_tests)
