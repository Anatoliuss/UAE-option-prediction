#include <catch2/catch_test_macros.hpp>
#include "black_scholes.hpp"
#include "monte_carlo.hpp"
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>  
using Catch::Approx;                 

TEST_CASE("MC converges to Black - Scholes (euro call, q=0)") {
    const double S=100, K=100, r=0.03, q=0.0, sigma=0.2, T=1.0;
    const double bs = call_price(S,K,r,sigma,T);

    auto [p1,se1] = mc_euro_option_price(S,K,r,q,sigma,T,true,  20000, 123, true);
    auto [p2,se2] = mc_euro_option_price(S,K,r,q,sigma,T,true, 100000, 123, true);

    // Within 3 standard errors is a common statistical acceptance band
    CHECK( std::fabs(p1 - bs) < 3.0 * se1 );
    CHECK( std::fabs(p2 - bs) < 3.0 * se2 );
}

TEST_CASE("MC put vs BS (q=0)"){
    const double S=120, K=100, r=0.01, q=0.0, sigma=0.25, T=0.5;
    const double bs = put_price(S,K,r,sigma,T);
    auto [p,se] = mc_euro_option_price(S,K,r,q,sigma,T,false, 120000, 777, true);
    CHECK( std::fabs(p - bs) < 3.0 * se );
}

TEST_CASE("MC handles edge cases T=0 and sigma=0") {
    // T=0 -> intrinsic
    {
        auto [p,se] = mc_euro_option_price(100,100,0.05,0.0,0.2,0.0,true, 10);
        CHECK( p == Approx(0.0).margin(1e-12) );
        CHECK( se == Approx(0.0).margin(1e-12) );
    }
    // sigma=0 -> forward, discounted intrinsic
    {
        auto [p,se] = mc_euro_option_price(100,100,0.05,0.0,0.0,1.0,true, 10);
        // deterministic: ST = 100*e^{0.05}; price = e^{-0.05} * max(ST-K,0)
        const double ST = 100*std::exp(0.05);
        const double disc = std::exp(-0.05);
        const double ref = disc * std::max(0.0, ST - 100.0);
        CHECK( p == Approx(ref).margin(1e-12) );
        CHECK( se == Approx(0.0).margin(1e-12) );
    }
}
