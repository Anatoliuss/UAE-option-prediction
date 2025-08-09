cmake_minimum_required(VERSION 3.25)
project(UAEOptionPricer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---- Dependencies: Catch2 ----
include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.4
)
FetchContent_MakeAvailable(Catch2)

# ---- Your code ----
# Header-only BS lives in include/, CSV loader is compiled
add_library(core STATIC
  src/csv_loader.cpp
)
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# (Optional) CLI app if you want to run src/main.cpp
add_executable(app src/main.cpp)
target_link_libraries(app PRIVATE core)

# ---- Tests ----
enable_testing()
add_executable(tests
  tests/test_bs.cpp
  tests/test_csv.cpp
  tests/test_greeks.cpp
  tests/test_vega.cpp
  tests/test_iv.cpp 
)
target_link_libraries(tests PRIVATE core Catch2::Catch2WithMain)
target_compile_definitions(tests PRIVATE PROJECT_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

include(CTest)
include(Catch)
catch_discover_tests(tests)
