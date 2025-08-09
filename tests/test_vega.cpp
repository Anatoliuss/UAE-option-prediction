#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "black_scholes.hpp"

TEST_CASE("Vega no-dividends") {
    // Case 1
    {
        const double S=100, K=100, r=0.05, sigma=0.2, T=1.0;
        // Expected ≈ 37.52403469169379
        REQUIRE(vega(S,K,r,sigma,T) == Catch::Approx(37.52403469169379).epsilon(1e-12));
    }
    // Case 2
    {
        const double S=100, K=110, r=0.05, sigma=0.2, T=1.0;
        // Expected ≈ 39.57604803881934
        REQUIRE(vega(S,K,r,sigma,T) == Catch::Approx(39.57604803881934).epsilon(1e-12));
    }
}