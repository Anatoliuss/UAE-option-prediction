#include <catch2/catch_test_macros.hpp>
#include "black_scholes.hpp"
#include <catch_test_macros.hpp>
#include <catch_amalgamated.hpp>

TEST_CASE("Call price textbook example") {
    double price = call_price(100.0, 100.0, 0.05, 0.20, 1.0);
    REQUIRE(price == Approx(10.4506).epsilon(1e-4));
}
