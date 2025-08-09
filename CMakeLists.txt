cmake_minimum_required(VERSION 3.25)
project(UAEOptionPricer LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.4          
)
FetchContent_MakeAvailable(catch2)

enable_testing()
add_executable(pricer-cli src/main.cpp)
target_include_directories(pricer-cli PRIVATE include)

# new tiny library
add_library(csv_loader src/csv_loader.cpp)
target_include_directories(csv_loader PUBLIC include)

# ---------------- Tests ----------------
add_executable(tests
    tests/test_csv.cpp
    tests/test_bs.cpp
)

target_link_libraries(tests
    PRIVATE
        csv_loader            # our CSV helper lib
        Catch2::Catch2WithMain
)

add_test(NAME all_unit_tests COMMAND tests)
