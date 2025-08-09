#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "CsvLoader.hpp"
#include <filesystem>
#include <string>

TEST_CASE("CSV loader reads two rows") {
    const std::string path = std::string(PROJECT_SOURCE_DIR) + "/data/data_raw/duae_sample.csv";
    REQUIRE(std::filesystem::exists(path));  // sanity check

    auto rows = load_csv(path);
    REQUIRE(rows.size() == 2);
    REQUIRE(rows[0].date == "2025-07-29");
    REQUIRE(rows[0].close == Catch::Approx(5114.99).epsilon(1e-9));
}
