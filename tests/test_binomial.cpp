#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "black_scholes.hpp"
#include "binomial.hpp"

TEST_CASE("Binomial converges to Black-Scholes (euro call, q=0)") {
    const double S=100, K=100, r=0.05, sigma=0.2, T=1.0;
    const double bs = call_price(S,K,r,sigma,T);
    for (int N : {10, 50, 100, 500}) {
        double b = binomial_price(S,K,r,sigma,T,N,true,false,0.0);
        REQUIRE(b == Catch::Approx(bs).epsilon(5e-2)); 
    }
}

TEST_CASE("American call with q=0 equals European") {
    const double S=100, K=100, r=0.05, sigma=0.2, T=1.0;
    for (int N : {25, 100, 300}) {
        double euro = binomial_price(S,K,r,sigma,T,N,true,false,0.0);
        double amer = binomial_price(S,K,r,sigma,T,N,true,true ,0.0);
        REQUIRE(amer == Catch::Approx(euro).epsilon(1e-12));
    }
}
TEST_CASE("American call with dividends >= European") {
    const double S=100, K=100, r=0.05, sigma=0.2, T=1.0, q=0.04;
    for (int N : {50, 200}) {
        double euro = binomial_price(S,K,r,sigma,T,N,true,false,q);
        double amer = binomial_price(S,K,r,sigma,T,N,true,true ,q);
        REQUIRE(amer >= euro);
    }
}
