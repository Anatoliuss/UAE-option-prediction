#include <catch2/catch_test_macros.hpp>
#include "CsvLoader.hpp"

TEST_CASE("CSV loader reads two rows") {
    auto rows = load_csv("data/data_raw/duae.csv");
    REQUIRE(rows.size() == 2);
    REQUIRE(rows[0].date == "2025-07-29");
    REQUIRE(rows[0].close == Approx(5114.99));
}
