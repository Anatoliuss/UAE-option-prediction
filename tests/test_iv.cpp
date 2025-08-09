#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "black_scholes.hpp"

TEST_CASE("IV (call) Newton+bisection recovers sigma across moneyness") {
    struct Case { double S,K,r,T,sigma; };
    const Case cases[] = {
        {100,100,0.05,1.0,0.2},   // ATM
        {100,120,0.05,1.0,0.35},  // OTM, higher vol
        {100, 80,0.05,2.0,0.15},  // ITM, longer T
    };
    for (auto c : cases) {
        const double mkt = call_price(c.S,c.K,c.r,c.sigma,c.T);
        const double iv  = implied_vol_call(c.S,c.K,c.r,c.T,mkt, /*sigma0=*/0.05);
        REQUIRE(iv == Catch::Approx(c.sigma).epsilon(1e-6));
    }
}
TEST_CASE("Implied vol (put) recovers true sigma") {
    const double S=100, K=100, r=0.05, T=1.0, true_sigma=0.25;
    const double market_price = put_price(S,K,r,true_sigma,T);
    const double iv = implied_vol_put(S,K,r,T, market_price, /*sigma0=*/0.1);
    REQUIRE(iv == Catch::Approx(true_sigma).epsilon(1e-6));
}
