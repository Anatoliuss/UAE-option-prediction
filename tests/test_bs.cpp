#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>  
#include "black_scholes.hpp"
#include <cmath>

TEST_CASE("Black-Scholes prices (no dividends)") { 
    const double S=100, K=100, r=0.05, sigma=0.2, T=1.0;
    const double c = call_price(S, K, r, sigma, T);
    const double p = put_price (S, K, r, sigma, T);
    REQUIRE(c == Catch::Approx(10.450583572185565).epsilon(1e-9));
    REQUIRE(p == Catch::Approx( 5.573526022256971).epsilon(1e-9));
    REQUIRE((c - p) == Catch::Approx(S - K*std::exp(-r*T)).epsilon(1e-10));
}

TEST_CASE("OTM/ITM (no dividends)") {
    const double S=100, K=110, r=0.05, sigma=0.2, T=1.0;
    const double c = call_price(S, K, r, sigma, T);
    const double p = put_price (S, K, r, sigma, T);
    REQUIRE(c == Catch::Approx(6.040088129724225).epsilon(1e-9));
    REQUIRE(p == Catch::Approx(10.675324824802779).epsilon(1e-9));
}